# SrO defect models
Thermodynamic defect models to extrapolate amount of SrO formed at operating conditions, from SrO formation energies calculated at T=0K by DFT.

## Dry air conditions
### Different cases

|  Case #      |  LSCF model  |  Origin of O  |  Origin of Sr  | Kröger-Vink defect equation                                                                             |
|:-----------: | :-----------:|:-----------:  | :-----------: |:-----------:                                                                                             |
| 1            | Slab        | O<sub>2(g)</sub>  | Surf       |$` Sr^{'}_{La} + \frac{1}{2} O_{2(g)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + 2B^{\bullet}_{B} `$       |
| 2            | Slab        | O<sub>2(g)</sub>  | Bulk       | same as 1                                                                                                |    
| 3            | Slab        | Bulk              | Surf       |$`Sr^{'}_{La} + O_O^X = SrO_s + V^{'''}_{La} + V_{O}^{\bullet \bullet}`$                                  |
| 4            | Slab        | Bulk              | Bulk       | same as 3 |
| 5            | Bulk        | O<sub>2(g)</sub>  | Bulk       |same as 1 |
| 6            | Bulk        | Bulk              | Bulk       |same as 3 |


### Defect model types

$`\Delta _{r} G(T,p) = 0 = \Delta _r G^{ref}(T,p) + RT ln (K)`$ for case 1 and 3.
$` \Delta _r G^{ref}(T,p)`$ is calcualted from formation energies (`dft_energies_OK.py` file in `./lib/`), by essentially replacing energies of gases by their respective chemical potentials at a given temperature $`T`$. 
For case 1

$`
\begin{align*}
K = \frac{\left[V^{'''}_{La}\right] \cdot \left[ B^{\bullet}_{B} \right ]^{2} }{\left[Sr^{'}_{La}\right] \cdot \left[B^{X}_{B}\right]^{2} \cdot \sqrt{\frac{p_{O_2}}{p}}} $= exp\left(\frac{-\Delta _r  G^{ref}}{RT}\right)
\end{align*}
`$

  
## Humid air conditions
### Different cases

|  Case #      |  LSCF model  |  Origin of O                  |  Origin of Sr  | Kröger-Vink defect equation                                                                             |
|:-----------: | :-----------:|:-----------:                  | :-----------:  |:-----------:                                                                                            |
| 1            | Slab        | H<sub>2</sub>O<sub>(g)</sub>   | Surf           | $` Sr^{'}_{La} + H_{2}O_{(g)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + H_{2(g)} + 2B^{\bullet}_{B} `$     |
| 2            | Slab        | H<sub>2</sub>O<sub>(g)</sub>   |Surf            | $` Sr^{'}_{La} + H_{2}O_{(g)} = SrO_{(s)} + (2H^+)^{'}_{La} `$                                             |
| 3            | Slab        | H<sub>2</sub>O<sub>(g)</sub>   | Surf           |  $` Sr^{'}_{La} +H_{2}O_{(g)} + B^{X}_{B} = SrO_{(s)} + (H^+)^{''}_{La} +  \frac{1}{2} H_{2(g)} + B^{\bullet}_{B} `$ |
| 4         | Slab        | H<sub>2</sub>O<sub>(g)</sub>   | Bulk           |  same as 1 |
### Defect model types


## Effect of hydroxylated surfaces
### Different cases

|  Case #      |  LSCF model  |  Origin of O                  |  Origin of Sr  | Kröger-Vink defect equation                                                                             |
|:-----------: | :-----------:|:-----------:                  | :-----------:  |:-----------:                                                                                            |
| 1            | Slab        | H<sub>2</sub>O<sub>(ads)</sub>   | Surf           | $` Sr^{'}_{La} + H_{2}O_{(ads)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + H_{2(g)} + 1s+ 2B^{\bullet}_{B} `$     |
| 2            | Slab        | H<sub>2</sub>O<sub>(ads)</sub>   |Surf            | $` Sr^{'}_{La} + H_{2}O_{(ads)} = SrO_{(s)} + (2H^+)^{'}_{La} + 1s `$                                             |
| 3            | Slab        | H<sub>2</sub>O<sub>(ads)</sub>   | Surf           |  $` Sr^{'}_{La} + H_{2}O_{(ads)} + B^{X}_{B} = SrO_{(s)} + (H^+)^{''}_{La} +  \frac{1}{2} H_{2(g)} + 1s + B^{\bullet}_{B} `$ |
