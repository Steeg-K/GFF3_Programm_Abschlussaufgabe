# GFF3_Programm_Abschlussaufgabe
GFF3 Programm Abschlussaufgabe von Moritz Hobein und Katharina Steeg
```
Programmaufruf:
python gff3_Analyzer_KSMH.py <your.gff3> <Flags>
```
Im Zuge der Abschlussarbeit unseres Bioinformatik Moduls an der Uni, haben wir dieses Programm zur Analyse und Auswertung von GFF3-Dateien geschrieben.
Es zählt die verschiedenen vorhandenen Feature-Typen der Datei und gibt die prozentuale Anzahl an. Es werden außerdem die COG Kategorien extrahiert und ausgewertet. Eine verkürzte Zusammenfassung der Gene (mit Name, Feature-Typ, Start, Stop, Produkt und Strang) wird als CSV Datei ausgegeben, genauso wie die Nukleotid - und Aminosäuresequenzen der entsprechenden Gene. 
Nach beliebigen Termen kann auch über den Flag -m <deineSuche> gesucht und deren Anzahl ausgegeben werden. Zudem ist es möglich, mit -e alle Zeilen, die den Suchterm enthalten, als Übersichts-CSV auszugeben. 
Über den Flag -o können auch Suchen definiert werden, die bestimmte Terme ausschließt. Z.B. kann mit -o CDS - nach allen CDS gesucht werden, die sich nicht auf dem "-"-Strang befinden. Die Angaben "+" und "-" beziehen sich dabei auf den jeweiligen Strang.   
Über die Zahlen 1-8 im -o Flag können die gff3-Spalten einzeln nach ihren Inhalten durchsucht werden. 


## Beispiel plot
![Codonusage_GCF_001050915_Graph](https://user-images.githubusercontent.com/80681798/111205837-f6c4a280-85c7-11eb-842b-5b57011bc8ba.png)
