# From Tague Lab

rhessys_diff = function(a, b) {

  bd = b$bd-a$bd
  bdg = b$bdg-a$bdg
  bdg.wy = b$bdg.wy-a$bdg.wy
  bdg.wyd = b$bdg.wyd-a$bdg.wyd
  bd.wy = b$bd.wy-a$bd.wy
  bd.wyd = b$bd.wyd-a$bd.wyd

  hd = b$hd-a$hd
  hdg = b$hdg-a$hdg
  hd.wy = b$hd.wy-a$hd.wy
  hd.wyd = b$hd.wyd-a$hd.wyd
  hdg.wy = b$hdg.wy-a$hdg.wy
  hdg.wyd = b$hdg.wyd-a$hdg.wyd

  zd = b$zd-a$zd
  zdg = b$zdg-a$zdg
  zd.wy = b$zd.wy-a$zd.wy
  zd.wyd = b$zd.wyd-a$zd.wyd
  zdg.wy = b$zdg.wy-a$zdg.wy
  zdg.wyd = b$zdg.wyd-a$zdg.wyd

  pd = b$pd-a$pd
  pdg = b$pdg-a$pdg
  pd.wy = b$pd.wy-a$pd.wy
  pd.wyd = b$pd.wyd-a$pd.wyd
  pdg.wy = b$pdg.wy-a$pdg.wy
  pdg.wyd = b$pdg.wyd-a$pdg.wyd

  cd = b$cd-a$cd
  cdg = b$cdg-a$cdg
  cd.wy = b$cd.wy-a$cd.wy
  cd.wyd = b$cd.wyd-a$cd.wyd
  cdg.wy = b$cdg.wy-a$cdg.wy
  cdg.wyd = b$cdg.wyd-a$cdg.wyd
  res=list(
    bd=bd, bdg=bdg, bd.wy=bd.wy, bdg.wy=bdg.wy, bd.wyd=bd.wyd, bdg.wyd=bdg.wyd,
    hd=hd, hdg=hdg, hd.wy=hd.wy, hdg.wy=hdg.wy, hd.wyd=hd.wyd, hdg.wyd=hdg.wyd,
    zd=zd, zdg=zdg, zd.wy=zd.wy, zdg.wy=zdg.wy, zd.wyd=zd.wyd, zdg.wyd=zdg.wyd,
    pd=pd, pdg=pdg, pd.wy=pd.wy, pdg.wy=pdg.wy, pd.wyd=pd.wyd, pdg.wyd=pdg.wyd,
    cd=cd, cdg=cdg, cd.wy=cd.wy, cdg.wy=cdg.wy, cd.wyd=cd.wyd, cdg.wyd=cdg.wyd)

  return(res)
}

aggregate.clean = function (x, use, func)
{
  ncolx = ncol(x)
  nvar = length(use)
  tmp = aggregate(x, by = use, func)
  tmp = tmp[, c((nvar + 1):(ncolx + nvar), seq(from = 1, to = nvar))]
  tmp
}
