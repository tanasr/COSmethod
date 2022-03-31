# COSmethod
Fourier-Cosinus methode (COS) für die Bewertung von Optionspreisen. BSc Thesis @ZHAW Winterthur (2017)

## Management summary
Die bedeutende Schwierigkeit in der Quantitative Finance ist die explizite Berechnung von Erwartungswerten, welche für die Preisfindung der Optionen benötigt werden. Was zu einfachen Formeln in der klassischen Vorgehensweise führt, wenn die zugrundeliegende Zufallsgrösse bspw. durch eine geometrische Brown'sche Bewegung modelliert wird, erweist sich als schwieriger in anspruchsvolleren Modellierungsansätzen.

Damit der Preis einer Option berechnet werden kann, wird die Wahrscheinlichkeitsdichtefunktion des zugrundeliegenden stochastischen Prozesses benötigt. Diese liegt für viele Modelle allerdings nicht in der analytischen Form vor. Die charakteristische Funktion hingegen existiert für fast alle relevanten Modelle, oft sogar in einfacher Form. Ist nun diese charakteristische Funktion bekannt, so kann aus ihr mit Hilfe der inversen Fourier-Transformation die Dichtefunktion gewonnen werden.

Das Ziel dieser Arbeit besteht darin, im Black-Scholes und im Heston Modell den Preis einer Europäischen Call Option und einer discretely monitored Barrier Option zu bewerten. Die Resultate werden anhand mehrerer Parameterstudien zusammengefasst und grafisch zu Illustrationszwecken dargestellt. Zusätzlich werden Aussagen betreffend der Konvergenz, Abweichung zur analytischen Lösung sowie Performance getroffen

Es wurde gezeigt, dass zwanzig Cosinusterme für den Algorithmus ausreichen, um den Optionspreis mit einer Abweichung zur analytischen Lösung von wenigen Zehntel zu erhalten. Dabei war die Rechenzeit vernachlässigbar klein. Ausserdem konnte anhand zweier unterschiedlich formulierten charakteristischen Funktionen der sogenannte Heston trap aufgezeigt werden. Die implementierte COS-Methode sollte mit der direkten Berechnung als auch indirekten über die Put-Call Parität anhand weiteren Parametern getestet und überprüft werden, ob die Instabilität verbessert werden kann. Dabei wäre die Auswirkung der Bestimmung der Integralgrenzen auf die Optionspreise interessant zu untersuchen.
