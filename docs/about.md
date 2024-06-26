# About
### How to cite
If you use ciclope please cite:
> Iori et al., (2023). Ciclope: micro Computed Tomography to Finite Elements. Journal of Open Source Software, 8(84), 4952, https://doi.org/10.21105/joss.04952 <br>

``` bibtex
@article{Iori2023, doi = {10.21105/joss.04952}, url = {https://doi.org/10.21105/joss.04952}, year = {2023}, publisher = {The Open Journal}, volume = {8}, number = {84}, pages = {4952}, author = {Gianluca Iori and Gianluigi Crimi and Enrico Schileo and Fulvia Taddei and Giulia Fraterrigo and Martino Pani}, title = {Ciclope: micro Computed Tomography to Finite Elements}, journal = {Journal of Open Source Software} }
```

### Statement of need
Micro Finite Element (microFE) models are often derived from volumetric stacks of micro Computed Tomography (microCT) images to non-destructively assess mechanical properties of biological or artificial specimens. However, the general absence of reproducible image-to-model pipelines and the use of proprietary or non-open-source software strongly limit validation and comparison of results across studies. **Ciclope** is a fully open-source pipeline, allowing to preprocess microCT data and obtain a corresponding microFE model, solve and postprocess it.

### Design
A typical pipeline for FE model generation from 3D microCT data is shown in the picture below.
![ciclope_design](ciclope_design.png)

* **Ciclope** is composed of a `core` module containing methods for voxel and tetrahedra FE model generation, and a module `utils` containing utilities for image and FE model pre- and post-processing. Both modules can be imported and used within Python. 
* Pipelines of FE model generation can be launched from the commandline using the `ciclope.py` script generated during installation. See the section [usage](usage) for more details.

### Ecosystem
**Ciclope** requires several dependencies. The following is a list of the main external packages required for FE model generation and solution:
* All mesh exports (voxel and tetrahedra Finite Elements) are performed with the [meshio](https://github.com/nschloe/meshio) module.
* Tetrahedra meshes are generated with [pygalmesh](https://github.com/nschloe/pygalmesh) (a Python frontend to [CGAL](https://www.cgal.org/)).
* High-resolution surface meshes for visualization are generated with the [PyMCubes](https://github.com/pmneila/PyMCubes) module.
* **Ciclope** generates Finite Element `.INP` files that can be solved using both [CalculiX](https://github.com/calculix) and [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/).
* The definition of material properties and of the FE analysis parameters (e.g. boundary conditions, simulation steps..) is handled through separate template files. The folders [material_properties](https://github.com/gianthk/ciclope/tree/master/material_properties) and [input_templates](https://github.com/gianthk/ciclope/tree/master/input_templates) contain a library of template files that can be used to generate FE simulations.
 * Additional libraries of [CalculiX](https://github.com/calculix) examples and template files can be found [here](https://github.com/calculix/examples) and [here](https://github.com/calculix/mkraska)
* For the post-processing of FE results, **Ciclope** uses [ParaView](https://www.paraview.org/) and the CalculiX to ParaView converter [`ccx2paraview`](https://github.com/calculix/ccx2paraview).