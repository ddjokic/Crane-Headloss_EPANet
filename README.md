# Crane-Headloss_EPANet
 
## Piping Hydraulic Calculation - Fluids based

Bassed on Caleb Bell's python library "fluids" - see manual. For fittings input refer to Fluids Documentation.

Calculates pressure drop in pipeline, based on Crane TP410m coefficients. Friction model is based on D'Arcy-Weisbach equation with Clemond's solution. On succesfull run scripts writes two text files in the folder (created by the script) named by pipe tag:

'results.txt' where all results of the run all collected, for flow step. Results include coefficients A, B, C in headloss equation (H[m]=AQ[cumh]**2+BQ[cumh}+c]);
EPAnet "crv" - file.
Q-H and Q-Re diagrams are ploted and located in the folder pipetag.

EPAnet Headloss "crv"-file for fittings only - to be used as input as headloss for EN's General Purpose Valve (GPV). In EN GPV dialog enter curve's EN-index as setting value. Note that Headloss in CRV-file is without friction, as EN calculates friction. For friction model in EN, use D-W (D'Arcy-Weisbach) equation for compatibility.

No warranty for the results!!!

Literature:

Clamond, Didier. “Efficient Resolution of the Colebrook Equation.” Industrial & Engineering Chemistry Research 48, no. 7 (April 1, 2009)

Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane, 2009

Idelchik, I. E, and A. S Ginevski ̆ı. Handbook of Hydraulic Resistance. Redding, CT: Begell House, 2007

Caleb Bell, Fluids Documentation, October 2017, http://pierreproulx.espaceweb.usherbrooke.ca/fluids.pdf, http://fluids.readthedocs.io/
