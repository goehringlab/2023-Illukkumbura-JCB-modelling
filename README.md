# Illukkumbura et al., 2023

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/goehringlab/2023-Illukkumbura-JCB-modelling/HEAD?filepath=%2Fscripts/notebook.ipynb)
[![CC BY 4.0][cc-by-shield]][cc-by]

Code for performing PDE modelling of advective transport in Illukkumbura et al., 2023

<p align="center">
    <img src="scripts/animation.gif" width="80%" height="80%"/>
</p>

## Instructions

See [this notebook](https://nbviewer.org/github/goehringlab/2023-Illukkumbura-JCB-modelling/blob/master/scripts/notebook.ipynb) for a demonstration of the model. To run the notebook interactively you have two options:

####  Option 1: Binder

To run in the cloud using Binder, click here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/goehringlab/2023-Illukkumbura-JCB-modelling/HEAD?filepath=%2Fscripts/notebook.ipynb)

(Please note that it may take several minutes to open the notebook)


#### Option 2: Docker

Step 1: Open [Docker](https://www.docker.com/products/docker-desktop/) and pull the docker image (copy and paste into the terminal)

    docker pull tsmbland/2023-illukkumbura-jcb-modelling

Step 2: Run the docker container (copy and paste into the terminal)

    docker run -p 8888:8888 tsmbland/2023-illukkumbura-jcb-modelling

This will print a URL for you to copy and paste into your web browser to open up Jupyter

Step 3: When finished, delete the container and image
    
    docker container prune -f
    docker image rm tsmbland/2023-illukkumbura-jcb-modelling


## Citation

This work is featured in the following article:

Rukshala Illukkumbura, Nisha Hirani, Joana Borrego-Pinto, Tom Bland, KangBo Ng, Lars Hubatsch, Jessica McQuade, Robert G. Endres, Nathan W. Goehring; Design principles for selective polarization of PAR proteins by cortical flows; J Cell Biol (2023) 222 (8): e202209111; doi: [https://doi.org/10.1101/2022.09.05.506621](https://doi.org/10.1083/jcb.202209111)

## License

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
