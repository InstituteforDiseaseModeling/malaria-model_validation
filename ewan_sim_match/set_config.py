def set_config( config, tmp_loc = [], rate= 1.0, infectivity = 1.0 ):
    config.parameters.Simulation_Type = "MALARIA_SIM" 
    # config.parameters.Acquisition_Blocking_Immunity_Decay_Rate = 0
    # config.parameters.Acquisition_Blocking_Immunity_Duration_Before_Decay = 0
    config.parameters.Incubation_Period_Exponential = 233.33
    config.parameters.Infectious_Period_Constant = 0
    config.parameters.Enable_Birth = 1
    #config.parameters.Enable_Coinfection = 1
    config.parameters.Enable_Demographics_Birth = 1
    config.parameters.Enable_Demographics_Reporting = 0
    # config.parameters.Enable_Immune_Decay = 0
    config.parameters.Migration_Model = "NO_MIGRATION"
    config.parameters.Run_Number = 99 
    config.parameters.Simulation_Duration = 70*365
    config.parameters.Enable_Demographics_Risk = 1
    config.parameters.Enable_Maternal_Infection_Transmission = 0
    config.parameters.Enable_Natural_Mortality = 1

    return config
