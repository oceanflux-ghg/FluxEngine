#Contains a dummy implementation of a custom gas transfer parameterisation used
# to test importing custom parameterisations from specified files

from fluxengine.core.rate_parameterisation import KCalculationBase



class example_custom_k_parameter_without_args(KCalculationBase):
    def __init__(self):
        pass
    
    def nput_names(self):
        return []
    
    def output_names(self):
        return []
    
    def __call__(self, data):
        return True


class example_custom_k_parameter_with_args(KCalculationBase):
    def __init__(self, customConfigVar):
        self.customConfigVar = customConfigVar
    
    def nput_names(self):
        return []
    
    def output_names(self):
        return []
    
    def __call__(self, data):
        return True