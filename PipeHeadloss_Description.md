##Piping Hydraulic Calculation - Fluids based

Bassed on <b>Caleb Bell's python library "fluids"</b> - see manual.
For fittings input refer to Fluids Documentation.

Calculates pressure drop in pipeline, based on <b>Crane TP410m</b> coefficients. Friction model is based on <b>D'Arcy-Weisbach</b> equation with Clemond's solution. On succesfull run scripts writes two text files in the folder (created by the script) named by pipe tag:
1. 'results.txt' where all results of the run all collected, for flow step. Results include coefficients A, B, C in headloss equation (H[m]=A*Q[cumh]**2+B*Q[cumh}+c]);
2. EPAnet "<b>crv</b>" - file.

Q-H and Q-Re diagrams are ploted and located in the folder <b>pipetag</b>.


<b>EPAnet</b> Headloss "<b>crv</b>"-file for fittings only - to be used as input as headloss for EN's General Purpose Valve (GPV). 
In EN GPV dialog enter curve's EN-index as setting value. Note that Headloss in CRV-file is <b>without friction</b>, as EN calculates friction.
For friction model in EN, use D-W (D'Arcy-Weisbach) equation for compatibility.

No warranty for the results!!!

<i>Literature:</i>

Clamond, Didier. “Efficient Resolution of the Colebrook Equation.” Industrial & Engineering Chemistry Research 48, no. 7 (April 1, 2009) 

Crane Co. Flow of Fluids Through Valves, Fittings, and Pipe. Crane, 2009 

Idelchik, I. E, and A. S Ginevski ̆ı. Handbook of Hydraulic Resistance. Redding, CT: Begell House, 2007 

Caleb Bell, Fluids Documentation, October 2017, http://pierreproulx.espaceweb.usherbrooke.ca/fluids.pdf, http://fluids.readthedocs.io/