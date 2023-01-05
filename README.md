# bulk-RNAseq-deconvolution-with-scRNAseq
The deconvolution for bulk RNAseq is to have a cell type fraction ratio, say differentational vs un-differential cell fraction in the sample.

Key packages :
-- DWLS
-- DeconRNASeq

The required input is the expression matrix of bulk RNAseq, and the single cell expression matrix with cell labels(characters with name of the cell types,Seurat RenameIdents function can assign the cell types based on the PanglaoDB/Cellmarkers database).

1) build cell signature file by the funcion buildSignatureMatrixMAST DWLS R package. This will be genes x cell types matrix, with the value of the percentages of the cell types.

2) Deconvolution bulk RNAseq by the function from DeconRNASeq R package. Output is sample x cell matrix, with the value of the percentages of the cell types.
