import math

# Alveolar Partial Pressure - Used in the Haldane Equation
def AlveolarPress(Q, PAmb0):
    """ Calculates the alveolar partial pressure of an inert gas at a given 
    ambient pressure. Uses a Respiratory Quotient of 0.7, which is the most
    conservative value. May be changed to the range 0.7-1.0. Inputs are Q =
    Inert Gas Quotient (FN2) and Ambient pressure at the beginning of
    the segment in bar. """

    RQ = 0.7      # Respiratory Quotient
    PH2O = 0.0627 # Water Vapor Pressure
    PCO2 = 0.0534 # CO2 Pressure in lungs

    return Q * (PAmb0 - PH2O + (1-RQ)/RQ * PCO2)

# Calculate Ceilings
def CalcCeiling(CompartmentLoadings, MValues):
    """ Calculates the ascent ceiling for a given compartment state and set
    of M-Values. """
    CompartmentCeilings = [0]*16
    for i in range(len(CompartmentLoadings)):
        CompartmentCeilings[i] = (CompartmentLoadings[i] - MValues[i][2]) / MValues[i][3]
    Ceiling = max(CompartmentCeilings)
    if Ceiling <= 0:
        return 0
    else:
        return Ceiling
		
# Static Haldane Equation
def SegmentStatic(CompartmentLoadings, PAmb0, Q, t, k):
    """ Calculates the inert gas loading of each tissue compartment for a
    given segment of the dive plan with a constant depth. Inputs are the 
    current Compartment Loadings in bar, ambient pressure at the beginning
    of the segment in bar, Inert Gas Quotient, duration of the segment in 
    mins, and the array of compartment k values. """

    Palv0 = AlveolarPress(Q, PAmb0)

    for i in range(len(CompartmentLoadings)):
        CompartmentLoadings[i] = Palv0 + (CompartmentLoadings[i] - Palv0) * math.exp(-k[i]*t)

    return CompartmentLoadings

# Dynamic Haldane Equation
def SegmentDynamic(CompartmentLoadings, PAmb0, Q, t, k, R):
    """ Calculates the inert gas loading of each tissue compartment for a given
    segment of the dive plan with linearly varying depth. Inputs are the 
    current Compartment Loadings in bar, ambient pressure at the beginning
    of the segment in bar, Inert Gas Quotient, duration of the segment in 
    mins, the array of compartment k values, and the rate of change in pressure,
    measured in bar/min. """

    Palv0 = AlveolarPress(Q, PAmb0)
    for i in range(len(CompartmentLoadings)):
        CompartmentLoadings[i] = Palv0 + R * (t - 1/k[i]) - (Palv0 - CompartmentLoadings[i] - R / k[i]) * math.exp(-k[i]*t)
    
    return CompartmentLoadings
