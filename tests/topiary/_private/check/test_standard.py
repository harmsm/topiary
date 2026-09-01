
import topiary
from topiary._private.check.standard import check_bool
from topiary._private.check.standard import check_float
from topiary._private.check.standard import check_int
from topiary._private.check.standard import check_iter
from topiary._private.check.standard import column_to_bool

import pytest
import numpy as np
import pandas as pd

def test_check_bool():

    true_values = [True,1.0,1,np.ones(1,dtype=np.bool_)[0]]
    for t in true_values:
        assert check_bool(t)

    false_values = [False,0.0,0,np.zeros(1,dtype=np.bool_)[0]]
    for f in false_values:
        assert not check_bool(f)

    bad_value = [None,123,-1,bool,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float]
    for b in bad_value:
        with pytest.raises(ValueError):
            value = check_bool(b)


def test_check_float():

    value = check_float(1.0)
    assert value == 1.0

    bad_value = [None,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float,np.nan]
    for b in bad_value:
        with pytest.raises(ValueError):
            value = check_float(b)

    good_value = [-np.inf,np.inf,0,1,"1.0"]
    for g in good_value:
        value = check_float(g)

    with pytest.raises(ValueError):
        check_float(1.0,minimum_allowed=2.0)

    with pytest.raises(ValueError):
        check_float(1.0,minimum_allowed=1.0,minimum_inclusive=False)

    value = check_float(1.0,minimum_allowed=1.0,minimum_inclusive=True)
    assert value == 1

    with pytest.raises(ValueError):
        check_float(1.0,maximum_allowed=0.5)

    with pytest.raises(ValueError):
        check_float(1.0,maximum_allowed=1.0,maximum_inclusive=False)

    value = check_float(1.0,minimum_allowed=1.0,maximum_inclusive=True)
    assert value == 1


def test_check_int():

    value = check_int(1)
    assert value == 1

    bad_value = [None,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float,int,np.inf,np.nan,1.3]
    for b in bad_value:
        print(b)
        with pytest.raises(ValueError):
            value = check_int(b)

    good_value = [-10,0,10,"10",10.0]
    for g in good_value:
        value = check_int(g)

    with pytest.raises(ValueError):
        check_int(1,minimum_allowed=2.0)

    with pytest.raises(ValueError):
        check_int(1,minimum_allowed=1,minimum_inclusive=False)

    value = check_int(1,minimum_allowed=1,minimum_inclusive=True)
    assert value == 1

    with pytest.raises(ValueError):
        check_int(1,maximum_allowed=0)

    with pytest.raises(ValueError):
        check_int(1,maximum_allowed=1,maximum_inclusive=False)

    value = check_int(1,minimum_allowed=1,maximum_inclusive=True)
    assert value == 1

_ITER_DF = pd.DataFrame({"test":[1,2,3]})

# Values check_iter must reject as not-an-iterable
BAD_ITERABLES = [None,0,list,float,int,np.inf,np.nan,1.3]

# Values check_iter must accept
GOOD_ITERABLES = [[],(),"test",{},_ITER_DF,np.arange(10)]


@pytest.mark.parametrize("bad",BAD_ITERABLES)
def test_check_iter_rejects_non_iterables(bad):

    with pytest.raises(ValueError):
        check_iter(bad)


@pytest.mark.parametrize("good",GOOD_ITERABLES)
def test_check_iter_accepts_iterables(good):

    # Note: this case used to be dead. The original test built a list of good
    # values and then looped over the *bad* list a second time, so nothing ever
    # checked that a valid iterable is accepted.
    out = check_iter(good)

    # check_iter hands the value back, so the caller can use its return
    assert out is not None or len(good) == 0


@pytest.mark.parametrize("value,required_type",[
    ([1,2,3],list),
    (tuple([1,2,3]),tuple),
    (np.arange(5),type(np.arange(1))),
    (_ITER_DF,type(_ITER_DF)),
])
def test_check_iter_accepts_matching_required_iter_type(value,required_type):

    out = check_iter(value,required_iter_type=required_type)
    assert isinstance(out,required_type)


@pytest.mark.parametrize("value,wrong_type",[
    ([1,2,3],str),   ([1,2,3],tuple),   ([1,2,3],type(_ITER_DF)),
    (tuple([1,2,3]),str), (tuple([1,2,3]),list), (tuple([1,2,3]),type(_ITER_DF)),
    (np.arange(5),str), (np.arange(5),tuple), (np.arange(5),list),
    (_ITER_DF,str), (_ITER_DF,tuple), (_ITER_DF,list),
])
def test_check_iter_rejects_wrong_required_iter_type(value,wrong_type):

    with pytest.raises(ValueError):
        check_iter(value,required_iter_type=wrong_type)


@pytest.mark.parametrize("value",[["test"],[1.0],[None],[1,1.0]])
def test_check_iter_rejects_wrong_required_value_type(value):

    with pytest.raises(ValueError):
        check_iter(value,required_value_type=int)


def test_check_iter_accepts_matching_required_value_type():

    assert np.array_equal(check_iter([1,2,3],required_value_type=int),[1,2,3])


@pytest.mark.parametrize("kwargs",[
    {"minimum_allowed":2},
    {"minimum_allowed":1,"minimum_inclusive":False},
])
def test_check_iter_rejects_values_below_minimum(kwargs):

    with pytest.raises(ValueError):
        check_iter([1],**kwargs)


@pytest.mark.parametrize("kwargs",[
    {"maximum_allowed":1},
    {"maximum_allowed":2,"maximum_inclusive":False},
])
def test_check_iter_rejects_values_above_maximum(kwargs):

    with pytest.raises(ValueError):
        check_iter([1,2],**kwargs)


def test_check_iter_accepts_inclusive_bounds():

    assert np.array_equal(check_iter([1],minimum_allowed=1,minimum_inclusive=True),[1])
    assert np.array_equal(check_iter([1,2],maximum_allowed=2,maximum_inclusive=True),[1,2])


@pytest.mark.parametrize("value,is_not_type",[
    ("test",str),
    ("test",[str,list]),
    ([1,2],[str,list]),
])
def test_check_iter_rejects_is_not_type(value,is_not_type):

    # In the original these three lived inside a single `with pytest.raises`
    # block, so only the first ever executed -- the other two were dead.
    with pytest.raises(ValueError):
        check_iter(value,is_not_type=is_not_type)


@pytest.mark.parametrize("value,is_not_type",[
    ([1],str),
    (tuple([1,2]),[str,list]),
])
def test_check_iter_accepts_allowed_is_not_type(value,is_not_type):

    assert np.array_equal(check_iter(value,is_not_type=is_not_type),value)


def test_check_iter_returns_the_value():

    assert np.array_equal(check_iter([1]),[1])


def test_column_to_bool():

    df = pd.DataFrame({"test":[True,False,
                               "True","False",
                               1,0,
                               "1","0",
                               "yes","no",
                               "T","F",
                               "Y","N"]})

    expected = np.array((1,0,1,0,1,0,1,0,1,0,1,0,1,0),dtype=bool)
    out = column_to_bool(df["test"],"test")
    assert np.array_equal(out,expected)

    df = pd.DataFrame({"test":[True,False,
                               "True","False",
                               1,0,
                               "1","0",
                               "yes","no",
                               "T","F",
                               "X","N"]})
    with pytest.raises(ValueError):
        out = column_to_bool(df["test"],"test")
