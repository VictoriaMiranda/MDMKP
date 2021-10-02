# MDMKP

Las instancias y su respectivas eficiencias se encuentran en las carpetas benchmark/250 y eficiencias, respectivamente.

El archivo media-worst_A-TSTS.cpp contiene el código de esta estrategia. 

Los parámetros que toma el programa son:

name: nombre de la instancia

time: tiempo de una ejecución (60 [s])

core: tamaño del core (probé con 1, 5, y 10, no probaría más allá de 20)

Itt: número de iteraciones sin mejora para calcular semi-opuesto [fase 2] (se multiplica por 1000, he probado con 12 y 15, probaría con valores un poco mas bajos como 10 o hasta 8) 

r = cantidad a restar en matriz tabú (he probado con 1 y 2, sería bueno probar hasta 5 al menos)

pct_L: tamaño de lista tabú, según porcentaje de N (antes no lo tenía como parámetro, sino que estaba fijo en 5 -> L = 0.05*N = 12)

lambda: factor que acompaña penalización en función de evaluación (también estaba fijo en 100)

seed: semilla que usará el programa.

Para ejecutar el código abro una cmd y me voy al path donde esten las carpetas benchmark, eficiencias, y media-worst_A-TSTS.cpp. Una vez en este path escribo "g++ -o main media-worst_A-TSTS.cpp". Luego ingreso "main name time core Itt r pct_L lambda seed". Por ejemplo, si los archivos estuviesen en una carpeta que se llame "A-TS", la cmd se vería algo así

C:\A-TS>g++ -o main media-worst_A-TSTSc.cpp

C:\A-TS>main 250-5-5-0-0.txt 60 1 12 1 5 100 0
