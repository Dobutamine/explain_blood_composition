from _blood_composition.lib import GetBloodComposition;


approx = GetBloodComposition(6.9, 27.1, 35.9, 25.0, 1.64, 0.0, 8.0, 5.0, 37.0);
print(approx.ph)
print(approx.pco2)
print(approx.hco3)
print(approx.be)
print(approx.po2)
print(approx.so2)
print(approx.steps_ab)
print(approx.steps_o2)