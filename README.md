# NestSurvivalSimulation

Contains several files that simulate multinomial nest fate data.

#"Sims_for_Tool.R" 
This is an R script that simulates nest fate data, analyzes the data using BUGS code, and then uses the resulting parameter estimates in a stochastic projection model to predict population growth rates with and without exclosures. The BUGS model and projection model are the same as those used in the tool, and the purpose of the script was to simulate data under varying sample sizes in order to determine the minimum sample size for tool use. There is a single "reference" run of 1000 nests, which is used to establish the correct management decision given the parameters used to simulate data. 
