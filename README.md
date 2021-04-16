# LNG Regasification

Python toolbox to model model the potential of LNG regasification heat recovery, in which the latent heat of LNG is used as heat sink for the cycle. The considered systems include CO2 transcritical and supercritical cycles, as well as Organic Rankine Cycle (ORC). The model is designed to converge in all conditions (subcritical, supercritical) and with all working fluids. The heat exchangers are modeled through a constant pinch point temperature difference.

Getting Started
---------------

Install CoolProp:

```bash
pip install CoolProp
```

Install Pyswarm:

```bash
pip install pyswarm
```

Run one of the simulation files, e.g. "Transcritical power cycle.py"


Example result
--------------
![R1336mzz(Z)](/docs/figures/fig.png)


Main authors
--------------
- Haoshui Yu (DTU, Denmark)
- Sylvain Quoilin (University of Liège, Belgium)

Reference
--------------

This toolbox was presented in the following publication:

Haoshui Yu1, Ning Guo, Sylvain Quoilin, Gürkan Sin, Performance comparison of organic Rankine cycle (ORC) and CO2 cycle for simultaneous utilization of LNG cold energy and solar energy, in Proceedings of the ORC2021 Conference, 2021.
