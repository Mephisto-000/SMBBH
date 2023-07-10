# SMBBH
The Study on the Velocity Ratio of Binary Black Holes with Galactic Potential.

In my thesis, we primarily study the relationship between velocity ratio and mass ratio in the binary black holes with galactic potential. In the field of astronomy, it is widely accepted that the center of active galactic nuclei hosts supermassive black hole. In the study of dual active galactic nuclei (AGNs), [Wang et al.](https://iopscience.iop.org/article/10.1088/0004-637X/705/1/L76) discovered a significant correlation between the velocity ratio and mass ratio of dual AGNs. Therefore, considering the influence of galactic potential, we explore the relationship between velocity ratio and mass ratio by choosing different initial conditions in binary black holes models.

## Platform of OS

**Ubuntu 22.04 LTS**

## Installation

Prerequisites:

- Anaconda 23.5.0
- Python 3.8

Install all dependencies in separate conda environment:
```shell
$ conda env create -f separate environment.yml
```

Activate this new environment:
```shell
$ conda activate smbbh
```

## Contents

$E_{0}$ : initial total energy.  
$V_{i}(0)$ : initial velocity of particle 1 or particle 2.  
$|v_{1,\hat{z}}|/|v_{2,\hat{z}}|$ : orbital radial velocity ratio.  
$m_{2}/m_{1}$ : mass ratio.  

### Results:  

<br/><br/>

| Model |  $E_{0}$   | $V_{i}^{2}(0)$ | $m_{2}/m_{1}$ | maximum $|v_{1,\hat{z}}|/|v_{2, \hat{z}}|$ | minimum $|v_{1,\hat{z}}|/|v_{2, \hat{z}}|$ |
| :---: | :--------: | :------------: | :-----------: | :----------------------------------------: | :----------------------------------------: |
|  A1   | $-25.1250$ |   $12.5625$    |     $1.0$     |                  $1.0000$                  |                  $1.0000$                  |
|  A2   | $-2.7917$  |    $1.3958$    |     $1.0$     |                  $1.0000$                  |                  $1.0000$                  |
|  B1   | $-25.1250$ |   $12.5625$    |     $1.5$     |         $1.5+6.3283\cdot 10^{-14}$         |         $1.5-7.6161\cdot 10^{-14}$         |
|  B2   | $-2.7917$  |    $1.3958$    |     $1.5$     |         $1.5+7.8826\cdot 10^{-14}$         |         $1.5-1.7675\cdot 10^{-13}$         |
|  C1   | $-25.1250$ |   $12.5625$    |   $2.3333$    |       $2.3333+2.0295\cdot 10^{-13}$        |       $2.3333-8.8374\cdot 10^{-14}$        |
|  C2   | $-2.7917$  |    $1.3958$    |   $2.3333$    |       $2.3333+1.5676\cdot 10^{-13}$        |       $2.3333-7.0166\cdot 10^{-14}$        |
|  D1   | $-25.1250$ |   $12.5625$    |     $4.0$     |         $4.0+4.4942\cdot 10^{-13}$         |         $4.0-6.0885\cdot 10^{-13}$         |
|  D2   | $-2.7917$  |    $1.3958$    |     $4.0$     |         $4.0+2.8244\cdot 10^{-13}$         |         $4.0-1.4033\cdot 10^{-13}$         |



### Model Demos:

- [Model A1](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_A1.ipynb)
- [Model A2](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_A2.ipynb)
- [Model B1](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_B1.ipynb)
- [Model B2](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_B2.ipynb)
- [Model C1](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_C1.ipynb)
- [Model C2](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_C2.ipynb)
- [Model D1](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_D1.ipynb)
- [Model D2](https://github.com/Mephisto-000/SMBBH/blob/main/demo/model_D2.ipynb)

## License

[MIT](https://github.com/Mephisto-000/SMBBH/blob/main/LICENSE)
