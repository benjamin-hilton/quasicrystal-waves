import scipy as sp

def generate_by_name(function, generations, mass_ratio, mass_default=1.0):
    return globals()[function](generations, mass_ratio, mass_default)

def fibonacci_1D(generations, mass_ratio, mass_default=1.0):
    current, previous = [mass_default], [mass_default * mass_ratio]

    for generation in xrange(generations):
        current, previous = current + previous, current

    return current

def thuemorse_1D(generations, mass_ratio, mass_default=1.0):
    current = sp.array([False])
    
    for generation in xrange(generations):
        current = sp.concatenate((current, ~current))

    current = current.astype(int) * mass_default*(mass_ratio - 1) #cast to floats
    current += mass_default

    return current

def equal_masses_1D(generations, mass_ratio, mass_default=1.0):
    return sp.ones(generations) * mass_default

def random_masses_1D(generations, mass_ratio, mass_default=1.0):
    random_list = [0.728399156869,0.938377471694,0.994319661004,0.586932514061,0.0385500165234,0.860967290702,0.844825329333,0.619469297834,0.47409338794,0.642620558302,0.762701405275,0.907921227632,0.360288502028,0.865239612806,0.933009037826,0.724257431181,0.115599875278,0.497313501387,0.486177616525,0.62977252722,0.413863185773,0.635731450165,0.376909079436,0.0452544551004,0.450749079233,0.686479056552,0.983443043149,0.42509664803,0.359036944948,0.642237551808,0.450672798637,0.00891825561225,0.438759261778,0.547109663094,0.0193541180154,0.0309856562098,0.842074421667,0.714447177067,0.987560796178,0.204514588001,0.457968658652,0.866465934442,0.792616307685,0.181426834936,0.996535435916,0.438135228257,0.977578639015,0.11669149113,0.705202812005,0.604835595026,0.78502049701,0.700794614141,0.834004684924,0.426081296434,0.135552785113,0.970141286627,0.929883042564,0.0829681152052,0.166707470083,0.32175102756,0.56308912989,0.378940402111,0.256787313697,0.0821909164531,0.609556484566,0.959328188558,0.545538674889,0.759279235013,0.208051604574,0.636085734117,0.0723232849794,0.190742689418,0.0953925758483,0.481400020394,0.90935206165,0.938668423622,0.944500176221,0.0595751411424,0.826466624481,0.808354938658,0.44793281137,0.150639536468,0.630222161601,0.894662844499,0.494848369408,0.28571471913,0.473982651151,0.0698983958875,0.515628632932,0.136261200407,0.965965673021,0.390182701682,0.178678232058,0.140925718493,0.605430162849,0.659958556773,0.235081815144,0.778360744217,0.953647414137,0.237915878167]
    return random_list[:generations]

def alternating_masses_1D(generations, mass_ratio, mass_default=1.0):
    current = sp.ones(generations) * mass_default
    current[::2] *= mass_ratio
    return current

def triangle_parity_1D(generations, mass_ratio, mass_default = 1.0):
    naturals = sp.arange(1, 1 + generations)
    triangles = (naturals) * (naturals + 1) * 0.5
    current = (triangles % 2) * mass_default * (mass_ratio - 1)
    current += mass_default
    return current

def arithmetic_1D(generations, mass_ratio, mass_default = 1.0):
    naturals = sp.arange(1, 1 + generations)
    current = []

    for index, natural in enumerate(naturals):
        even = index % 2 == 0

        if even:
            current += list(sp.ones(natural) * mass_default * mass_ratio)
        else:
            current += list(sp.ones(natural) * mass_default)

    return current
