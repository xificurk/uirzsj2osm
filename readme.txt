Skripty pro import sídelních jednotek z UIR-ZSJ do OSM.

Příprava před importem:
	Potřebné knihovny a programy:
		python 2.x: http://python.org/download/
		pyproj: http://code.google.com/p/pyproj/downloads/list
	Stáhnutí UIR-ZSJ:
		Do adresáře uir-zsj rozbal obsah archivu http://notes.czso.cz/csu/rso.nsf/i/uir111dc/$File/uir111dc.zip (nebo jeho novější verzi stáhnutou z http://www.czso.cz/csu/rso.nsf/i/prohlizec_uir_zsj). Jediné potřebné soubory jsou OBCE.DBF, COBE.DBF a ZSJD.DBF.
	V adresáři se skripty vytvořte soubor credentials.txt, do kterého na první řádek napište své uživatelské jméno a na druhý řádek heslo pro přihlášení ke svému OSM účtu.

Přesný postupu při importu je na wiki: http://wiki.openstreetmap.org/wiki/Import_UIR-ZSJ
