The main plotting files are in root directory with appropriate names.
The different defect models are in library folders with the suffix models in it.

## Defect models for LSCF decomposition
Thermodynamic defect models to extrapolate amount of SrO formed at operating conditions, from SrO formation energies calculated at T=0K by DFT.

## SrO formation

### Dry air conditions

|  Case #      |  LSCF model  |  Origin of O      |  Origin of Sr  | Kröger-Vink defect equation                                                                              |
|:-----------: | :-----------:|:-----------:      | :-----------:  |:-----------:                                                                                             |
| 1            | Slab         | O<sub>2(g)</sub>  | Surf           |$` Sr^{'}_{La} + \frac{1}{2} O_{2(g)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + 2B^{\bullet}_{B} `$       |
| 2            | Slab         | O<sub>2(g)</sub>  | Bulk           | same as 1                                                                                                |    
| 3            | Slab         | Bulk              | Surf           |$`Sr^{'}_{La} + O_O^X = SrO_s + V^{'''}_{La} + V_{O}^{\bullet \bullet}`$                                  |
| 4            | Slab         | Bulk              | Bulk           | same as 3                                                                                                |
| 5            | Bulk         | O<sub>2(g)</sub>  | Bulk           |same as 1                                                                                                 |
| 6            | Bulk         | Bulk              | Bulk           |same as 3                                                                                                 |


  
### Humid air conditions

|  Case #      |  LSCF model  |  Origin of O                  |  Origin of Sr  | Kröger-Vink defect equation                                                                                          |
|:-----------: | :-----------:|:-----------:                  | :-----------:  |:-----------:                                                                                                         |
| 1            | Slab         | H<sub>2</sub>O<sub>(g)</sub>  | Surf           | $` Sr^{'}_{La} + H_{2}O_{(g)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + H_{2(g)} + 2B^{\bullet}_{B} `$               |
| 2            | Slab         | H<sub>2</sub>O<sub>(g)</sub>  | Surf           | $` Sr^{'}_{La} + H_{2}O_{(g)} = SrO_{(s)} + (2H^+)^{'}_{La} `$                                                       | 
| 3            | Slab         | H<sub>2</sub>O<sub>(g)</sub>  | Surf           |  $` Sr^{'}_{La} +H_{2}O_{(g)} + B^{X}_{B} = SrO_{(s)} + (H^+)^{''}_{La} +  \frac{1}{2} H_{2(g)} + B^{\bullet}_{B} `$ |
| 4            | Slab         | H<sub>2</sub>O<sub>(g)</sub>  | Bulk           |  same as 1                                                                                                           |


### Hydroxylated surfaces

|  Case #      |  LSCF model  |  Origin of O                  |  Origin of Sr  | Kröger-Vink defect equation                                                                                                     |
|:-----------: | :-----------:|:-----------:                  | :-----------:  |:-----------:                                                                                                                    |
| 1            | Slab         | H<sub>2</sub>O<sub>(ads)</sub>   | Surf           | $` Sr^{'}_{La} + H_{2}O_{(ads)} + 2B^{X}_{B} = SrO_{(s)} + V^{'''}_{La} + H_{2(g)} + 1s+ 2B^{\bullet}_{B} `$                 |
| 2            | Slab         | H<sub>2</sub>O<sub>(ads)</sub>   |Surf            | $` Sr^{'}_{La} + H_{2}O_{(ads)} = SrO_{(s)} + (2H^+)^{'}_{La} + 1s `$                                                        |
| 3            | Slab         | H<sub>2</sub>O<sub>(ads)</sub>   | Surf           |  $` Sr^{'}_{La} + H_{2}O_{(ads)} + B^{X}_{B} = SrO_{(s)} + (H^+)^{''}_{La} +  \frac{1}{2} H_{2(g)} + 1s + B^{\bullet}_{B} `$ |


## Volatile Sr(OH)2 formation


|  Case #      |  LSCF model  |  Origin of O            |  Origin of Sr  | Kröger-Vink defect equation                                                                                                     |
|:-----------: | :-----------:|:-----------:            | :-----------:  |:-----------:                                                                                                                    |
| 1            | Slab         | Slab and atmosphere     | Surf           | $` Sr^{'}_{La} + H_{2}O_{(ads)} + O^{X}_{O} = Sr(OH)_{2(g)} + V^{'''}_{La} +  V^{\bullet \bullet}_{O}  `$                 |
