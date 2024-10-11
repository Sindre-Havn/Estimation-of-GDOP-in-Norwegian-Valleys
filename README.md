# Estimation-of-GDOP-in-Norwegian-Valleys
This project was my summerjobb of 2024 in the norwegian road-department (ITS avdelingen i Statens vegvesen).

The project was to devolope a procedure to calculate GDOP along road where nearby moutains shield field of view to gnss satellites, specifically the norwegian valley Romsdalen. This repo serves as protype solution to that task.

This is done by using "Skyline" tool in ArcGIS Pro and public available GNSS ephemeris data.
This solution has a know bug where about half of the possible skylines fail to be generated, reason unkown.

The pdf: "Beregning_av_DOP_langs_norsk_vei_med_fjellskygge_effekter" is my notebook troughout the project.

The pdf: "Sommerjobb vegvesenet 2024" is the final project presentation.



To get started type these in terminal:

Create virtual environment:             python -m venv .venv

Activate environment windows:           .venv\Scripts\activate.bat.
                for other OS:           .venv/bin/activate.


Install packages with:                  pip install -r requirements.txt