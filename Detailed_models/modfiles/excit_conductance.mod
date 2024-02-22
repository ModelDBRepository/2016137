
TITLE Channel: excit_conductance

COMMENT
    Simple example of a leak/passive conductance. Note: for GENESIS cells with a single leak conductance,
        it is better to use the Rm and Em variables for a passive current.
ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX excit
    ? A non specific current is present
    RANGE e
    NONSPECIFIC_CURRENT i
    
    RANGE gmax, gion, i_excit, syn_per_area, tau, rate
    
}

PARAMETER { 
      

    gmax = 0.8 (nS)  ? default value, should be overwritten when conductance placed on cell
    
    e = 0 (mV) ? default value, should be overwritten when conductance placed on cell
    
    tau = 0.0127 (s)
    
    rate = 3 (Hz)
    
    syn_per_area = 1 (1/um2) ? default value, should be read from a separate dat file in the main model
    
}



ASSIGNED {
      

    v (mV)
        
    i (mA/cm2)
    i_excit (mA/cm2)
    gion (S/cm2)
        
}

BREAKPOINT { 

    gion = gmax * tau * rate * syn_per_area * 0.1 ? multiply by 0.1 because of unit conversion
    i = gion * (v - e)
    i_excit = i
        

}


