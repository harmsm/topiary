
import pytest

import topiary

from topiary.raxml.model import _generate_parsimony_tree
from topiary.raxml.model import _model_thread_function
from topiary.raxml.model import _parse_raxml_info_for_aic
from topiary.raxml.model import find_best_model

import pandas as pd

import os, json

def test__parse_raxml_info_for_aic(tmpdir,small_phylo):
    """
    _parse_raxml_info_for_aic is a pure parser over a raxml log file, so test it
    against both a real log written by raxml and hand-built edge cases.
    """

    # 1. A real raxml log from the committed test data
    real_log = small_phylo["01_gene-tree/working/infer-ml-tree/alignment.phy.raxml.log"]
    out = _parse_raxml_info_for_aic(real_log)

    assert out["L"] == pytest.approx(-294.266461)
    assert out["N"] == 21
    assert out["AIC"] == pytest.approx(630.532923)
    assert out["AICc"] == pytest.approx(1092.532923)
    assert out["BIC"] == pytest.approx(655.272053)

    # 2. Hand-built minimal file -- same keys, different values
    info_file = os.path.join(tmpdir,"info.log")
    with open(info_file,"w") as f:
        f.write("some preamble we do not care about\n")
        f.write("Final LogLikelihood: -100.5\n")
        f.write("AIC score: 201.0 / AICc score: 202.0 / BIC score: 203.0\n")
        f.write("Free parameters (model + branch lengths): 7\n")

    out = _parse_raxml_info_for_aic(info_file)
    assert out == {"L":pytest.approx(-100.5),
                   "AIC":pytest.approx(201.0),
                   "AICc":pytest.approx(202.0),
                   "BIC":pytest.approx(203.0),
                   "N":7}

    # 3. A log with none of the interesting lines yields an empty dict rather
    #    than raising
    empty_file = os.path.join(tmpdir,"empty.log")
    with open(empty_file,"w") as f:
        f.write("nothing to see here\n")

    assert _parse_raxml_info_for_aic(empty_file) == {}

    # 4. A missing file is an error
    with pytest.raises(FileNotFoundError):
        _parse_raxml_info_for_aic(os.path.join(tmpdir,"not-a-file.log"))


def test__generate_parsimony_tree_builds_expected_command(mocker,tmpdir):
    """
    _generate_parsimony_tree is a thin wrapper that assembles a raxml command,
    so check the command rather than running raxml.
    """

    run_raxml = mocker.patch("topiary.raxml.model.run_raxml")

    _generate_parsimony_tree("alignment.phy",
                             run_directory="parsimony-tree",
                             seed=True,
                             num_threads=1,
                             raxml_binary="raxml-ng")

    run_raxml.assert_called_once()
    kwargs = run_raxml.call_args.kwargs

    assert kwargs["algorithm"] == "--start"
    assert kwargs["alignment_file"] == "alignment.phy"
    assert kwargs["run_directory"] == "parsimony-tree"
    assert kwargs["num_threads"] == 1
    assert kwargs["raxml_binary"] == "raxml-ng"
    assert "--tree pars{1}" in " ".join(kwargs["other_args"])


def test__model_thread_function_returns_parsed_scores(mocker,tmpdir):
    """
    _model_thread_function runs raxml for one model and then parses the log it
    produced. Mock the raxml call and hand it a log to parse.
    """

    os.chdir(tmpdir)

    def _fake_run_raxml(**kwargs):
        run_dir = kwargs["run_directory"]
        os.makedirs(run_dir,exist_ok=True)
        with open(os.path.join(run_dir,"alignment.phy.raxml.log"),"w") as f:
            f.write("Final LogLikelihood: -50.0\n")
            f.write("AIC score: 110.0 / AICc score: 120.0 / BIC score: 130.0\n")
            f.write("Free parameters (model + branch lengths): 5\n")

    mocker.patch("topiary.raxml.model.run_raxml",side_effect=_fake_run_raxml)

    kwargs = {"alignment_file":"alignment.phy",
              "tree_file":"tree.newick",
              "model":"LG",
              "run_directory":"model-LG",
              "seed":12345,
              "num_threads":1,
              "raxml_binary":"raxml-ng"}

    out = _model_thread_function(kwargs)

    assert out["L"] == pytest.approx(-50.0)
    assert out["N"] == 5
    assert out["AIC"] == pytest.approx(110.0)
    assert out["AICc"] == pytest.approx(120.0)
    assert out["BIC"] == pytest.approx(130.0)

    # The run directory is cleaned up on the way out
    assert not os.path.exists("model-LG")


def test__model_thread_function_returns_none_when_raxml_crashes(mocker,tmpdir):
    """
    A raxml crash for one model must not kill the whole model search -- the
    thread returns None and find_best_model skips that model.
    """

    os.chdir(tmpdir)

    def _crash(**kwargs):
        os.makedirs(kwargs["run_directory"],exist_ok=True)
        raise RuntimeError("raxml fell over")

    mocker.patch("topiary.raxml.model.run_raxml",side_effect=_crash)

    kwargs = {"alignment_file":"alignment.phy",
              "tree_file":"tree.newick",
              "model":"LG",
              "run_directory":"model-LG",
              "seed":12345,
              "num_threads":1,
              "raxml_binary":"raxml-ng"}

    assert _model_thread_function(kwargs) is None


@pytest.mark.run_raxml
def test_find_best_model(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]
    os.chdir(tmpdir)

    model_matrices = ["LG","JTT"]
    model_rates = ["","G8"]
    model_freqs = ["","FO"]
    model_invariant = ["","IO"]
    calc_dir = "test_out"

    find_best_model(df,
                    model_matrices=model_matrices,
                    model_rates=model_rates,
                    model_freqs=model_freqs,
                    model_invariant=model_invariant,
                    calc_dir=calc_dir)

    # Make sure we're building models properly
    all_models = []
    for a in model_matrices:
        for b in model_rates:
            for c in model_freqs:
                for d in model_invariant:
                    model = [m for m in [a,b,c,d] if m != ""]
                    all_models.append("+".join(model))

    # Read output dataframe
    out_df = pd.read_csv(os.path.join("test_out","output","model-comparison.csv"))

    # Make sure it is the right length, sorted correctly, and has all models
    assert len(out_df) == len(all_models)
    assert out_df["p"].iloc[0] > out_df["p"].iloc[-1]
    assert len(set(all_models).intersection(set(out_df["model"]))) == len(all_models)

    # Make sure we can read sane json
    f = open(os.path.join("test_out","run_parameters.json"))
    out_json = json.load(f)
    f.close()

    # Make sure best model is recorded properly in output json
    best_model = out_df["model"].iloc[0]
    assert out_json["model"] == best_model
    assert out_json["calc_status"] == "complete"

