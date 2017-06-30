# evolving-tradeoffs
Code and data to reproduce the results of "Identifying and modelling 
evolving tradeoffs in multipurpose water resources systems via agent-based
negotiation protocols".

# Usage
The code works in MATLAB R2017a
* First navigate to the base folder of the repository, and run INSTALL.m
* You might want to check the system definition, and the stepcosts. Details
  on the system equations, and on how the simulation works, are contained in 
  code/Test/Test_simulateSystem.m. Parameters can be found in 
  code/Test/Test_environment.m and code/Test/Test_parameters.m;
* You might also want to generate a synthetic scenario. Run 
  code/script_decisionScenarioBuilding.m, it will call the function to 
  generate the inflows (code/Test/Test_generatePhysicalScenario.m) and uses
  the policies identified via EMODPS to compose a decision scenario, such as
  data/decisionScenario.mat;
* To run negotiations, you can choose to exhaustively explore the attitudes
  of the agents: in this case you need code/script_fewNegotiations.m;
* Otherwise, you might want to adopt the availability bias-based attitude 
  model and run code/script_negotiation.m to exhaustively explore the effect 
  of the memory parameter, the only parameter of the attitude model;
* The equation of the attitude model can be found in the class
  CooperativeAgent, function updateConcessionCoefficient() in the file 
  code/Negotiation/CooperativeAgent.m;

# Credits
Thanks to Paolo Gazzotti for his first version of the code. Thanks to 
article coauthors for advices and ideas, Matteo Giuliani (mxgiuliani00), 
Andrea Castelletti and Francesco Amigoni. This is their work as well.