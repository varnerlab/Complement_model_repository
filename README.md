## Reduced order complement model

### Background ###
Complement is an important pathway in innate immunity, inflammation, and many disease processes.
However, despite its importance, there have been few validated mathematical models of complement activation.
In this study, we developed an ensemble of experimentally validated reduced order complement models.
We combined ordinary differential equations with logical rules to produce a compact yet predictive complement model.
The model,  which described the lectin and alternative pathways, was an order of magnitude smaller than comparable models in the literature.
We estimated an ensemble of model parameters from in vitro dynamic measurements of the C3a and C5a complement proteins.
Subsequently, we validated the model on unseen C3a and C5a measurements not used for model training.
Despite its small size, the model was surprisingly predictive.
Global sensitivity and robustness analysis suggested complement was robust to any single therapeutic intervention.
Only the knockdown of both C3 and C5 consistently reduced C3a and C5a formation from all pathways.
Taken together, we developed a reduced order complement model that was computationally inexpensive,
and could easily be incorporated into pre-existing or new pharmacokinetic models of immune system function.
The model described experimental data, and predicted the need for multiple points of therapeutic intervention to fully disrupt complement activation.

The complement model is described in the publication:

[Sagar A, Dai W, Minot M and J. Varner (2016) Reduced order modeling and analysis of the human complement system. bioRxiv 059386; doi: http://dx.doi.org/10.1101/059386](http://biorxiv.org/content/early/2016/06/16/059386)

### Installation
You can download this repository as a zip file, or clone or pull it by using the command:

	git pull https://github.com/varnerlab/Complement_model_repository

or

	git clone https://github.com/varnerlab/Complement_model_repository

### Model code and parameter ensemble
The complement model equations were implemented in [Julia](http://julialang.org) and solved using the CVODE routine of the [Sundials package](https://github.com/JuliaLang/Sundials.jl). The model code and parameter ensemble is freely available under an [MIT software license](https://opensource.org/licenses/MIT).

The model equations are encoded in ``Balances.jl`` which is called by the ``SolveBalances.jl`` driver function. The user should not directly call ``SolveBalances.jl``. Rather, multiple parameter sets can be simulated by calling the driver function from a script. The kinetic and other model parameters are encoded in ``DataFile.jl`` as a dictionary. The parameters stored in this dictionary can be updated in memory to run different simulations. An example script to simulate the model over the parameter ensemble is encoded in ``sample_ensemble.jl``. To execute this script, issue the command:

``julia> include("sample_ensemble.jl")``

This will load the parameter ensemble from disk (assumed to be stored in the ``data`` subdirectory), update the parameter dictionary returned by ``DataFile.jl`` with the new parameters, and solve the model equations. The Pareto rank of the simulation can be adjusted by changing the selection criteria (L13):

	# Select the desired rank -
	idx_rank = find(rank_array .<= 5.0)

This code will simulate all parameter sets in the ensemble with Pareto rank five or less.

The model ensemble was estimated using the [JuPOETs package](https://github.com/varnerlab/POETs.jl). The complement objective functions, and problem specific constraints are encoded in ``complement_lib.jl``. The JuPOETs package is described in the publication:

[Bassen D, Vilkhovoy M, Minot M, J. Butcher and J. Varner (2016) JuPOETs: A Constrained Multiobjective Optimization Approach to Estimate Biochemical Model Ensembles in the Julia Programming Language. bioRxiv 056044; doi: http://dx.doi.org/10.1101/056044](http://biorxiv.org/content/early/2016/05/30/056044)

__Prerequisites__: [Julia](http://julialang.org) and the [Sundials package](https://github.com/JuliaLang/Sundials.jl) must be installed on your computer before the model equations can be solved. In addition, in the example routine ``sample_ensemble.jl`` the ensemble output is plotted using the [PyPlot](https://github.com/stevengj/PyPlot.jl) package which requires a working Python installation.  

### Test simulation routines ###
To test your model and [Julia](http://julialang.org) installation we have included routines to recreate Fig 2 and 3 of the [Sagar et al study](http://biorxiv.org/content/early/2016/06/16/059386). In each of these routines, we sample the ensemble and plot the mean and 95% CI of
the ensemble versus the corresponding experimental data.

Filename | Figure | Readout | zymosan (mg/ml)
--- | --- | --- | ---
``Make_Fig_2A.jl`` | Fig. 2A | C3a | 0.0
``Make_Fig_2B.jl`` | Fig. 2B | C5a | 0.0
``Make_Fig_2C.jl`` | Fig. 2C | C3a | 1.0
``Make_Fig_2D.jl`` | Fig. 2D | C5a | 1.0
``Make_Fig_3_C1_C3a.jl`` | Fig. 3 | C3a | 0.1
``Make_Fig_3_C1_C5a.jl`` | Fig. 3 | C5a | 0.1
``Make_Fig_3_C2_C3a.jl`` | Fig. 3 | C3a | 0.01
``Make_Fig_3_C2_C5a.jl`` | Fig. 3 | C5a | 0.01
``Make_Fig_3_C3_C3a.jl`` | Fig. 3 | C3a | 0.001
``Make_Fig_3_C3_C5a.jl`` | Fig. 3 | C5a | 0.001

### Data directory ###
In this study we reproduced data from the study of [Shaw and coworkers](https://www.ncbi.nlm.nih.gov/pubmed/?term=Morad+and+Shaw+2015):

[Morad HO, Belete SC, Read T and AM Shaw (2015) Time-course analysis of C3a and C5a quantifies the coupling between the upper and terminal Complement pathways in vitro. J Immunol Methods. 2015 Dec;427:13-8. doi: 10.1016/j.jim.2015.09.001](https://www.ncbi.nlm.nih.gov/pubmed/?term=Morad+and+Shaw+2015)

to train and test the effective complement model. The data used for model training and validation is contained in the ``data`` subdirectory.

Filename | Original Figure | Current Figure | Species | Role
--- | --- | --- | --- | ---
``Shaw2015_Fig2a_C3a.txt`` | Fig. 2A | Fig. 2A | C3a | training
``Shaw2015_Fig3ai_C5a_original.txt`` | Fig. 3a(i) | Fig. 2B | C5a | training
``Shaw2015_Fig2e_C3a.txt`` | Fig. 2E | Fig. 2C | C3a | training
``Shaw2015_Fig3c_C5a_original.txt`` | Fig. 3C | Fig. 2D | C5a | training
``Shaw2015_Fig2d_C3a.txt`` | Fig. 2D | Fig. 3 (1,1) | C3a | prediction
``Shaw2015_Fig3b_C5a_original.txt`` | Fig. 3B	| Fig. 3 (2,1) | C5a | prediction
``Shaw2015_Fig2c_C3a.txt`` | Fig. 2C | Fig. 3 (1,2) | C3a | prediction
``Shaw2015_Fig3aiii_C5a_original.txt`` | Fig. 3A(iii) | Fig. 3 (2,2) | C5a | prediction
``Shaw2015_Fig2b_C3a.txt`` | Fig. 2B | Fig. 3 (1,3) | C3a | prediction
``Shaw2015_Fig3aii_C5a_original.txt`` | Fig. 3A(ii) | Fig. 3 (2,3) | C5a | prediction
