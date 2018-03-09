<img align="right" width="270" height="240"
     src="https://github.com/ruhugu/sandpiles/raw/master/output_examples/example_config.png">

# sandpiles
     
Implementation of the [Bak-Tang-Wiesenfeld (BTW) sandpiles model](https://en.wikipedia.org/wiki/Abelian_sandpile_model).
This was the first discovered dynamical system that showed criticality without the fine tuning of any of its parameters.
This property is called [self-organized criticality](https://en.wikipedia.org/wiki/Self-organized_criticality).

This python module allows to create, analyse and visualize lattices following this model.


## Required packages

This modules requires both [matplotlib](https://matplotlib.org/) and [numpy](http://www.numpy.org/) to be installed.


## Some output examples

Period 2 limit cycle in a 20x20 lattice with periodic boundary conditions ([script](sandpiles/scripts/limitcycledist.py)).

<img src="https://github.com/ruhugu/sandpiles/blob/master/output_examples/lcycle_random.png" alt="Drawing" width="500"/>

Cascade evolution ([script]( sandpiles/scripts/clusterevolution.py )):

<p class="indented"><img src="https://github.com/ruhugu/sandpiles/blob/master/output_examples/clusterevolutionL50.gif" alt="Drawing" width="350"/></p>

Cascade size and period (power law) distributions for several lattice sizes (scripts: [1](sandpiles/scripts/cascadestatistics.py) and [2](sandpiles/scripts/replotcascadestatistics.py)):

<img src="https://github.com/ruhugu/sandpiles/raw/master/output_examples/cascadeduration.png" alt="Drawing" width="600"/>
<img src="https://github.com/ruhugu/sandpiles/raw/master/output_examples/cascadesize.png" alt="Drawing" width="600"/>


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
