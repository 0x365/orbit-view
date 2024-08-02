# Orbit View

### About

A tool to visualise the stability of 3-body orbits.

The stability is given as standard distribution of the euclidean distance from the first orbit over all other orbits. The period of the first orbit is calculated by looking at the euclidean distance from the intial positions of the 3-bodys and using scipy.signal.find_peaks to identify close contacts with the intial positions.

Calculating the period and stability for all the 200x200 orbits in the grid is computationally expensive and is done with a different repository and is therefore not shown here currently.

This is just a tool to visualise and shouldn't be used as an official way to calculate orbit stability or orbit periods.

### Period Map and Stability Map
<img src="https://github.com/0x365/orbit-view/blob/main/data/period_map.png" width="300" height="250"></img>
<img src="https://github.com/0x365/orbit-view/blob/main/data/stability.png" width="300" height="250"></img>

### Setup
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Run
```bash
python main.py
```

Then just click on the plot to view the orbit at that specific stability
