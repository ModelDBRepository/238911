from Motility_Model import *


controlStim = ExtracellStim("Control")
control = MotilityModel_CES(controlStim,intensity="moderate")

control['Fig'].show()