# Code for simulating data using linear causal structure
# adapted in part from https://onlinelibrary.wiley.com/doi/full/10.1111/ele.14033
library(simstudy)
library(dagitty)
library(ggdag)
library(ggplot2)

#create DAG in dagity
dag1 <- dagitty('dag {
bb="-0.5,-0.5,0.5,0.5"
C1 [pos="0.049,-0.354"]
C2 [pos="-0.143,-0.405"]
C3 [pos="-0.343,-0.280"]
C4 [pos="-0.368,-0.045"]
C5 [pos="-0.128,0.102"]
E [exposure,pos="-0.205,-0.091"]
O [outcome,pos="0.085,-0.108"]
C1 -> O
C2 -> E
C2 -> O
C3 -> C1
C4 -> C3
C4 -> C5
C4 -> E
E -> C1
E -> C5
E -> O
O -> C5
}')

tidy_dagitty(dag1)
ggdag(dag1)
dagitty::adjustmentSets(dag1, effect = 'total')
dagitty::adjustmentSets(dag1, effect = 'direct')
ggdag_adjustment_set(dag1)

# set seed
set.seed(5) 
#using simstudy to simulate data using linear formulas in conjunction with the DAG

simdef <- defData(varname = "C4", dist = "normal", formula = 0, variance = 1)
simdef <- defData(simdef, varname = "C2", dist = "normal", formula = 0, variance = 1)
simdef <- defData(simdef, varname = "C3", dist = "normal", formula = "0.2 * C4", variance = 1)
simdef <- defData(simdef, varname = "E", dist = "normal", formula = "0.3 * C4 + 0.6 * C2", variance = 1)
simdef <- defData(simdef, varname = "C1", dist = "normal", formula = "-1.1 * E + -0.4 * C3", variance = 1)
simdef <- defData(simdef, varname = "O", dist = "normal", formula = "0.5 * C1 + -0.2 * E + 0.5 * C2 + -0.2 * C4", variance = 1)
simdef <- defData(simdef, varname = "C5", dist = "normal", formula = "0.5 * C4 + 0.5 * E + 0.5 * O", variance = 1)

# create 10000 observations 
simdata <- genData(10000,simdef)

# linear regression models of all variables vs adjustment set needed for causal inference 
noncausal_model<- glm(O ~ E + C1 + C2 + C3 + C4 + C5, data= simdata) 
summary(noncausal_model) # will return forestry estimate as well as AIC value 
BIC(noncausal_model) # will return BIC value 
confint(noncausal_model) # will give 95% CI for forestry estimate

#gives appropriate estimate for total effect of exposure on outcome
causal_model.1 <- glm(O ~ E + C2 + C4 , data= simdata) 
summary(causal_model.1) # will return forestry estimate as well as AIC value of model
BIC(causal_model.1)
confint(causal_model.1) # will give 95% CI for forestry estimate

#gives appropriate estimate for total effect of exposure on outcome
causal_model.2 <- glm(O ~ E + C2 + C3, data= simdata) 
summary(causal_model.2) # will return forestry estimate as well as AIC value of model
BIC(causal_model.2)
confint(causal_model.2) # will give 95% CI for forestry estimate

#direct effects model
causal_model.3 <- glm(O ~ E + C1 + C2 , data= simdata) 
summary(causal_model.3) # will return forestry estimate as well as AIC value of model
BIC(causal_model.3)
confint(causal_model.3) # will give 95% CI for forestry estimate

#NOTES
#effects are multipled between nodes in a path and 
#added to the values of converging paths to determine the total effects  

#additional model----
set.seed(5) 
#using simstudy to simulate data using linear formulas in conjunction with the DAG
simdef <- defData(varname = "E", dist = "normal", formula = 0, variance = 1)
simdef <- defData(simdef, varname = "C1", dist = "normal", formula = "2 * E", variance = 1)
simdef <- defData(simdef, varname = "C2", dist = "normal", formula = "7 * C1", variance = 1)
simdef <- defData(simdef, varname = "O", dist = "normal", formula = "2 * E + 5 * C2 + 3 * C1", variance = 1)

simdata <- genData(10000,simdef)

causal_model.4 <- glm(O ~ E + C1, data= simdata) 
summary(causal_model.4) # will return forestry estimate as well as AIC value of model
BIC(causal_model.4)
confint(causal_model.4) # will give 95% CI for forestry estimate



library(dagR)
confound_dag <- dag.init(
  # REQUIRED: SPECIFY COVARIATES AND THEIR TYPE
  covs = c(1,1,1,1,1) ,# specify that single covariate (common cause) is "known"
  
  # OPTIONAL: PROVIDE CUSTOM LABELS FOR DAG NODES AND SIM DATA
  symbols = c('E',  #0
              'O',  #1
              'C1', #2
              'C2', #3
              'C3', #4
              'C4', #5
              'C5'),#6 or -1
  
  # REQUIRED: SPECIFY ARCS BETWEEN NODES
  arcs = c(0, 1,   # arc from primary exposure to outcome
           0, 2,   # 
           2, 1,
           3, 0,
           3, 1,
           4, 2,
           5, 0,
           5 ,4,
           5, 6,
           0 ,6,
           1, 6,
           3 , 6),  # 
  
  # OPTIONAL: PROVIDE NAMES FOR X, C, AND O (USED IN LEGEND)
  x.name = 'forst_spA',   # name of primary exposure
  y.name = 'forst_spY',            # name of outcome
  cov.names = c('Common cause'), # (vector of) names of covariates
  
  # REQUIRED: TURN OFF WEIRD DEFAULTS
  noxy   = T        # tell dagR we want to specify our *own* arc between X and O
)
plot(confound_dag)

simdata_dagR <- dag.sim(
  # REQUIRED: PROVIDE DAG OBJECT AND DESIRED SAMPLE SIZE
  confound_dag,
  n = 10000,
  # REQUIRED: SPECIFY NODE TYPES (1 = BINARY, 0 = CTS)
  binary = c(0,  
             0,  
             0,
             0,
             0,
             0),
  
  # REQUIRED: SPECIFY AND MEAN AND SD OF CTS NODES
  # (*) MU = PREVALENCE FOR BINARY NODES
  mu = c(0,    # prevalence p of primary exposure X
         0,    # mean of common cause C
         0),    # "    "  outcome O
  stdev = c(1,   # sd of random noise added when generating X
            1,  # "  "  common cause C
            1),  # "  "  random noise added when generate O
  
  # REQUIRED: SPECIFY DIRECT EFFECTS OF ARCS
  b = c(-0.2,       # direct effect b of X -> O 
        -1.1,  # direct effect *exp(b)* of C -> X (since X binary)
        0.5),    # direct effect b of C -> O
  
  # OPTIONAL: SET SEED FOR REPRODUCIBILITY
  seed = 12345,  # set seed for reproducibility
  naming = 2,    # use custom labels for output variable names
  verbose = F    # suppress alg. output
)

Causal_Model3 <- glm(species_Y ~ forestry, data=confound_sim) 
summary(Causal_Model3) # will return forestry estimate as well as AIC value of model
BIC(Causal_Model3)
confint(Causal_Model3) # will give 95% CI for forestry estimate


#other resources
# https://iyarlin.github.io/2019/07/23/mixed_dag_simulation_using_simmixeddag_package/