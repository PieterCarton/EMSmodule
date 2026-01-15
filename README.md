# EMSmodule

<img src="images/FLXNconcept.png" alt="FLEXINet Conceptual Schematic" width="400"/>

**‘The aim of FLEXINet is a system that accelerates the energy transition. We hope to make a substantial contribution to reaching climate targets by cleverly combining various techniques – think of blending recycled batteries with flexible heat pumps and the charging of electric cars.’**

## Description
The code implements the Energy Management System (EMS), of Results 5 and 6 of the project. The EMS is a software that controls the energy flows in the building, in order to minimize the energy costs. The EMS is implemented in `Julia`, and it is based on [JuMP.jl](https://jump.dev/) and [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl). The EMS is an optimization based controller and it sends the optimal setpoints to the different devices.

The workflow is the following:

![EMS Workflow](images/grabs_AppEnergy.png)

## Badges
PENDING

## Installation
The common installation is followed. If necessary, add directly this repository to the Julia environment.

```julia
julia> ]
Pkg> add https://https://github.com/DarioSlaifsteinSk/EMSmodule.git # main
using EMSmodule
```

If you have any issues with the installation, please contact the developers.

## Usage

The functions that compose the EMS Library are in the files:
```julia
makeEMSobjs.jl # creates the mutable structs/objects
EMSfns.jl # functions to create the optimization model
ESSfns.jl # functions to add the ESS to the optimization models
EMSrfofns.jl # functions for RFO (Random Field Optimization)
thermfns.jl # thermal devices functions
makeEMSplots.jl # plotting functions
makeForecasts.jl # forecasting functions
simEMS.jl # simulation functions
testEMS.jl # testing functions UNDER DEVELOPMENT
```
Examples on how to use the package are under development.

## Sequential Decision Making

Following the [Universal Modelling Framework (UMF)](https://castle.princeton.edu/sda/) this library implements Direct Lookahead (DLA) Policies. The DLA is a model-based policy that uses a model of the system to predict the future and optimize the control actions. The available DLAs are a day-ahead (DA) planner and a Model Predictive Controller (MPC). The basic algorithm is depicted in the following figure:


![Approximated Direct Lookahead (DLA)](images/appDLA.png)


The DA planner uses a model of the system to predict the future and optimize the control actions for the next 24 hours. The MPC uses a model of the system to predict the future and optimize the control actions for the next 24 hours, but it also uses the actual measurements to update the model and the optimization problem every hour. The MPC is an economic non-linear MPC (NLP-eMPC) receding horizon controller.

## Forecasting

The forecasts can either be provided to the optimization model as input data, or they can be simulated by adding noise to the deterministic inputs.
A native `Julia` forecast package is under development.

## Support
For bugs or feature requests, please open an issue. For usage questions or extensions please refer to the discussion section.

<!-- ## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README. -->
<!-- 
## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser. -->

## Authors and acknowledgment

This repository contains the work produced for the [FLEXINET](https://www.tudelft.nl/en/eemcs/flexinet) project by the [DCE&S group](https://www.tudelft.nl/en/eemcs/the-faculty/departments/electrical-sustainable-energy/dc-systems-energy-conversion-storage), Electrical Sustainable Energy Dept. of the TU Delft. The work belongs to Dario Slaifstein, Gautam Rituraj, and Joel Alpizar.

## License
For open source projects, say how it is licensed. PENDING

<!-- ## Project status
Under development. -->

# References

This library is based and used in the following publications:


[1] D. Slaifstein, G. R. Chandra Mouli, L. Ramirez-Elizondo, and P. Bauer, “Sequential operation of residential energy hubs using physics-based economic nonlinear MPC,” Appl Energy, vol. 408, p. 127402, Apr. 2026, doi: 10.1016/J.APENERGY.2026.127402.

[2] D. Slaifstein, G. R. C. Mouli, L. Ramirez-Elizondo, and P. Bauer, “Aging-aware Energy Management for Residential Multi-Carrier Energy Systems,” Aug. 2025, Accessed: Mar. 22, 2025. [Online]. Available: http://arxiv.org/abs/2503.16139 *Under Review in Journal of Energy Storage*

[3] D. Slaifstein, A. M. Agudin, G. R. C. Mouli, L. Ramírez-Elizondo, and P. Bauer, “Stochastic Mobility Integration into Residential Energy Hubs,” in 2024 IEEE International Conference on Electrical Systems for Aircraft, Railway, Ship Propulsion and Road Vehicles &amp; International Transportation Electrification Conference (ESARS-ITEC), IEEE, Nov. 2024, pp. 1–7. doi: 10.1109/ESARS-ITEC60450.2024.10819794.

Cite as:
```bibtex
@article{Slaifstein2026,
   author = {Darío Slaifstein and Gautham Ram Chandra Mouli and Laura Ramirez-Elizondo and Pavol Bauer},
   doi = {10.1016/J.APENERGY.2026.127402},
   issn = {0306-2619},
   journal = {Applied Energy},
   month = {4},
   pages = {127402},
   publisher = {Elsevier},
   title = {Sequential operation of residential energy hubs using physics-based economic nonlinear MPC},
   volume = {408},
   url = {https://linkinghub.elsevier.com/retrieve/pii/S0306261926000541},
   year = {2026}
}
```

<!-- ## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.tudelft.nl/dces/emsmodule.git
git branch -M main
git push -uf origin main
```
## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.tudelft.nl/dces/emsmodule/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing (SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

## Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thanks to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README

Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.
-->
