import numpy as _np_
"""
.. module:: unbiasedSymmetricMetrics
   :platform: Unix
   :synopsis: Tool to calculate unbiased symmetric metrics.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

# library of unbiased symmetrics metrics.
#	Based on:"New unbiased symmetric metrics for evaluation of air quality models"
#	S. Yu et al.,Atmos. Sci. Let. 7: 26-34 (2006) DOI: 10.1002/asl.125
#	http://onlinelibrary.wiley.com/doi/10.1002/asl.125/pdf


def __calcSi__(m, o):
    """
	Calcuate S_i factor array. 

	"""
    si = lambda m, o: (m - o) / abs(m - o)
    try:
        s = (m - o) / _np_.ma.abs(m - o)
    except:
        s = []
        for mo, ob in zip(m, o):
            if mo == ob: s.append(1.)
            elif _np_.ma.masked in [mo, ob]: s.append(_np_.ma.masked)
            else:
                s.append(si(mo, ob))
        s = _np_.ma.array(s)
    return s


def __calcS__(m, o):
    """
	Calcuate S factor array. 
	"""
    return (m.mean() - o.mean()) / abs(m.mean() - o.mean())


def __testMNFBarrays__(model, obs, key):
    """ 
	Standard set of tests to ensure that the mean normalised factor bias/error calculates are stable.
	"""
    if len(model) != len(obs):
        print(key, "calculation:\tModel and data don't match in length!",
              len(model), '!=', len(obs))
        assert False
    if len(model) == len(obs) == 1:
        return _np_.ma.array(model), _np_.ma.array(obs)
    model = _np_.float64(_np_.ma.array(model))
    obs = _np_.float64(_np_.ma.array(obs))
    return model, obs


def MNFB(model, obs):
    """ 
	Calculate mean normalised factor bias. 
	Usage: MNFB(model, obs)
	"""
    model, obs = __testMNFBarrays__(model, obs, 'MNFB')
    si = __calcSi__(model, obs)
    return _np_.ma.sum(si *
                       (_np_.ma.exp(_np_.ma.abs(_np_.ma.log(model / obs))) -
                        _np_.ma.ones(len(obs)))) / float(len(obs))


def MNAFE(model, obs):
    """ 
	Calculate mean normalised absoulte factor error.
	Usage: MNAFE(model, obs)
	"""
    model, obs = __testMNFBarrays__(model, obs, 'MNAFE')
    return _np_.ma.sum(
        _np_.ma.abs(
            _np_.ma.exp(_np_.ma.abs(_np_.ma.log(model / obs))) -
            _np_.ma.ones(len(obs)))) / float(len(obs))


def NMBF(model, obs):
    """ 
	Calculate normalised mean bias factor.
	Usage: NMBF(model, obs)	
	"""
    model, obs = __testMNFBarrays__(model, obs, 'NMBF')
    s = __calcS__(model, obs)
    return s * (
        _np_.ma.exp(_np_.ma.abs(_np_.ma.log(model.sum() / obs.sum()))) - 1.)


def NMAEF(model, obs):
    """ 
	Calculate normalised mean bias factor.
	Usage: NMAEF(model, obs)	
	"""
    model, obs = __testMNFBarrays__(model, obs, 'NMAEF')
    s = __calcS__(model, obs)
    return _np_.ma.abs(model - obs).sum() / (obs.sum()**((1. + s) / 2.) *
                                             (model.sum()**((1. - s) / 2.)))
