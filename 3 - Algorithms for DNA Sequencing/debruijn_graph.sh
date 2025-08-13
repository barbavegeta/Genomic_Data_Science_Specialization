digraph "DeBruijn graph" {
 ATT [label="ATT"] ;
 AAG [label="AAG"] ;
 GAT [label="GAT"] ;
 TCT [label="TCT"] ;
 AGA [label="AGA"] ;
 TTC [label="TTC"] ;
 CTC [label="CTC"] ;
 CTA [label="CTA"] ;
 AAG -> AGA ;
 AGA -> GAT ;
 GAT -> ATT ;
 ATT -> TTC ;
 TTC -> TCT ;
 TCT -> CTC ;
 CTC -> TCT ;
 TCT -> CTA ;
}
