import scipy as sp

def pseudo_vectorised(output_types="", excluded=None, cache=True):
    """Vectorise the function using scipy.vectorize"""
    
    def _wrapper_func(func):
        return sp.vectorize(func, otypes=output_types, excluded=excluded, cache=cache)

    return _wrapper_func

