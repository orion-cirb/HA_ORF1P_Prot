* **Developed for:** Rania
* **Team:** Prochiantz
* **Date:** July 2023
* **Software:** Fiji

### Images description

3D images taken with a x25 objective.

3 channels:
  1. *405:* Nuclei (mandatory)
  2. *488:* Protein (mandatory)
  3. *642:* HA-ORF1p (optional)
     
If *.roi* file is provided, analysis is performed in it. Otherwise, analysis is performed in the entire image. 

### Plugin description

* Find best focused slice and detect nuclei in 2D with Cellpose
* If HA-ORF1p channel provided, find best focused slice and detect cells in 2D with Cellpose 
* Compute colocalization between nuclei and HA-ORF1p cells, tag each nucleus as HA-ORF1p+ or HA-ORF1p-
* For each cell, measure nucleus circularity + protein intensity in nucleus, inner nucleus, inner ring and outer ring

### Dependencies

* **Cellpose** conda environment + *cyto2* model
* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin

### Version history

Version 2 released on July 17, 2023.

