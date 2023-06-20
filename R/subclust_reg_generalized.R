# Note: gamma beta are actually transposed from the paper
# library(Rfast); library(mclust); library(Rcpp)

# Helper Functions --------------------------------------------------------

get.xyz = function(data, Gamma, px, py) {
  xyz = data %*% Gamma
  z_index = if (px + py == nrow(Gamma)) c() else (px + py + 1):nrow(Gamma)
  return(list(X = xyz[, 1:px, drop = FALSE],
              Y = xyz[, px + 1:py, drop = FALSE],
              Z = xyz[, z_index, drop = FALSE]))
}

rwts = function(ncol, nrow) {
  w = matrix(rexp(ncol * nrow), nrow, ncol)
  w = w / rowSums(w)
  return(w)
}

update.xyz = function(state) {
  xyz  = get.xyz(state$data, state$Gamma, state$px, state$py)
  state$X = xyz$X
  state$Y = xyz$Y
  state$Z = xyz$Z
  return(state)
}

update.bars = function(state) {
  z_ngh = state$z_ngh
  setx = state$setx
  sety = state$sety
  px = state$px
  py = state$py
  p = nrow(state$Gamma)

  state$mu.bar.g.h = lapply(1:state$G, function(k) {
    c(state$mu.g[, setx[k]],
      c(1, state$mu.g[, setx[k]]) %*% state$beta.g.h[, , sety[k]],
      rep(0, p - px - py))
  })

  # B.dim = dim(state$beta.g.h)[1:2] - c(1L, 0L)
  # state$sigma.bar.g.h.inv = lapply(1:state$G, function(k) {
  #   Sginv = solve(state$sigma.g[, , setx[k]])
  #   Sghinv = solve(state$sigma.g.h[, , sety[k]])
  #   B = state$beta.g.h[-1, , sety[k]]
  #   if (is.vector(B)) dim(B) = B.dim
  #   neg.BSghinv = -B %*% Sghinv
  #
  #   tmp = matrix(0, p, p)
  #   tmp[1:px, 1:px] = Sginv - tcrossprod(neg.BSghinv, B)
  #   tmp[1:px, px + 1:py] = neg.BSghinv
  #   tmp[px + 1:py, 1:px] = t(neg.BSghinv)
  #   tmp[px + 1:py, px + 1:py] = Sghinv
  #
  #   if (p > px + py) tmp[(px + py + 1):p, (px + py + 1):p] = diag(1/state$psi, p - px - py)
  #   return(tmp)
  # })

  state$sigma.bar.g.h.inv = asplit(
    update_bars_sigma_inv(state$G,
                          state$setx, state$sety,
                          state$sigma.g, state$sigma.g.h, state$beta.g.h,
                          state$psi,
                          px, py, p)
    , 3)
  return(state)
}

number.of.parameters = function(state) {
  reg = state$regression
  sub = state$subcluster
  px = state$px
  py = state$py
  Gx = state$Gx
  Gy = state$Gy
  p = nrow(state$Gamma)

  if (!is.null(state$diagonal) && state$diagonal) {
    npars.x = 2 * Gx * px
    npars.y = if (sub) {2 * sum(Gy) * py} else {2 * Gy[1] * py}
  } else {
    npars.x = (Gx - 1) * (px + px * (px + 1) / 2) + (px + px)
    npars.y = if (sub) {
      (sum(Gy) - 1) * (py + py * (py + 1) / 2) + (py + py)
    } else {
      (Gy[1] - 1) * (py + py * (py + 1) / 2) + (py + py)
    }
  }
  npars.beta = if (reg) {px * py * (if (sub) sum(Gy) else Gy[1])} else {0}
  npars.z = 2 * (p - px - py)
  npars.Gamma = p * (p - 1) / 2
  npars.tau = if (sub) state$G - 1 else Gx + Gy[1] - 2
  npars.total = npars.x + npars.y + npars.z + npars.Gamma + npars.beta + npars.tau

  return(npars.total)
}

get.bic = function(state) {
  if (length(state) == 1) return(NA)
  ll = tail(na.omit(state$ll), 1)
  return(number.of.parameters(state) * log(nrow(state$data)) - 2 * ll)
}

# Initialization ----------------------------------------------------------

init.mclust = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = svd(scale(data))$v
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  modelNames = if (diagonal) c("EII", "VII", "EEI", "VEI", "EVI", "VVI") else NULL
  if (subcluster) {
    z_ngh = Mclust(data, G, modelNames, verbose = FALSE)$z
  } else {
    zx = Mclust(xyz$X, Gx, modelNames, verbose = FALSE)$z
    zy = Mclust(xyz$Y, Gy[1], modelNames, verbose = FALSE)$z
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zy)
    }
  }

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster, diagonal = diagonal,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state, diagonal = diagonal)
  state$psi = state$psi * 0 + 10
  state = update.bars(state)

  for (i in 1:100) {
    state = e.step(state)
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }

  return(state)
}

init.mclust.randgamma = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = pracma::randortho(p)
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  modelNames = if (diagonal) c("EII", "VII", "EEI", "VEI", "EVI", "VVI") else NULL
  if (subcluster) {
    z_ngh = Mclust(data, G, modelNames, verbose = FALSE)$z
  } else {
    zx = Mclust(xyz$X, Gx, modelNames, verbose = FALSE)$z
    zy = Mclust(xyz$Y, Gy[1], modelNames, verbose = FALSE)$z
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zy)
    }
  }

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster, diagonal = diagonal,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state, diagonal = diagonal)
  state$psi = state$psi * 0 + 10
  state = update.bars(state)

  for (i in 1:100) {
    state = e.step(state)
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }

  return(state)
}

init.mclust.2 = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = TRUE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = svd(scale(data))$v
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  if (subcluster) {
    zx.model = Mclust(xyz$X, Gx, verbose = FALSE)
    zx = zx.model$z
    zys = lapply(1:Gx, function(g) {
      zy.model = Mclust(xyz$Y[zx.model$classification == g,], Gy[g], verbose = FALSE)
      return(predict(zy.model, xyz$Y)$z)
    })
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zys[[g]])
    }
  } else {
    zx = Mclust(xyz$X, Gx, verbose = FALSE)$z
    zy = Mclust(xyz$Y, Gy[1], verbose = FALSE)$z
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zy)
    }
  }

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = colMeans(xyz$Z^2),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state)
  state = update.bars(state)

  # for (i in 1:100) {
  #   state = e.step(state)
  #   state$Gamma = m_step_gamma_dual(
  #     state$G,
  #     state$sigma.bar.g.h.inv,
  #     state$mu.bar.g.h,
  #     state$z_ngh,
  #     state$data,
  #     state$Gamma
  #   )
  # }

  return(state)
}

init.mclust.3 = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = TRUE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = svd(data)$v
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  if (subcluster) {
    mc.X = Mclust(xyz$X, Gx, verbose = FALSE)
    mc.XYZ = Mclust(data %*% Gamma, G, verbose = FALSE)
    transpo = table(mc.X$classification, mc.XYZ$classification)
    perms = gtools::permutations(Gx, Gx)

    fit = lapply(1:nrow(perms), function(x) {
      lpSolve::lp.transport(
        transpo, "max",
        col.signs = rep("==", G), col.rhs = rep(1, Gx),
        row.signs = rep("==", Gx), row.rhs = Gy[perms[x,]]
      )
    })
    best.perm = perms[which.max(sapply(fit, `[[`, "objval")), ]

    zx = mc.X$z[,best.perm]
    perm.class = max.col(zx)
    zys = lapply(1:Gx, function(g) {
      zy.model = Mclust(xyz$Y[mc.X$classification == g,], Gy[g], verbose = FALSE)
      return(predict(zy.model, xyz$Y)$z)
    })
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zys[[g]])
    }
  } else {
    zx = Mclust(xyz$X, Gx, verbose = FALSE)$z
    zy = Mclust(xyz$Y, Gy[1], verbose = FALSE)$z
    z_ngh = NULL
    for (g in 1:Gx) {
      z_ngh = cbind(z_ngh, zx[,g] * zy)
    }
  }

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = colMeans(xyz$Z^2),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state)
  state = update.bars(state)

  # for (i in 1:100) {
  #   state = e.step(state)
  #   state$Gamma = m_step_gamma_dual(
  #     state$G,
  #     state$sigma.bar.g.h.inv,
  #     state$mu.bar.g.h,
  #     state$z_ngh,
  #     state$data,
  #     state$Gamma
  #   )
  # }

  return(state)
}

init.mclust.semi.sup = function(data, px, py, Gx, Gy, class = rep(NA, nrow(data)),
                                regression = TRUE, subcluster = TRUE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = svd(scale(data))$v
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  z_ngh = MclustSSC(data, class, G, verbose = FALSE)$z

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state)
  state$psi = state$psi * 0 + 10
  state = update.bars(state)

  z.held = state$z_ngh[!is.na(class),]
  for (i in 1:100) {
    state = e.step(state)
    state$z_ngh[!is.na(class),] = z.held
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }

  return(state)
}

init.mclust.semi.sup.z = function(data, px, py, Gx, Gy, class = rep(NA, nrow(data)),
                                regression = TRUE, subcluster = TRUE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = svd(scale(data))$v
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  z_ngh = MclustSSC(data, class, G, verbose = FALSE)$z

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state)
  state$psi = state$psi * 0 + 10
  state = update.bars(state)

  z.held = state$z_ngh[!is.na(class),]
  for (i in 1:100) {
    state = e.step(state)
    state$z_ngh[!is.na(class),] = z.held
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }

  return(state)
}

init.fancy = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE) {
  stopifnot(px >= 1)
  stopifnot(py >= 1)
  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)
  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  mc = Mclust(data, G, verbose = FALSE)
  Gamma = svd(t(mc$parameters$mean))$v
  rotmu = t(mc$parameters$mean) %*% Gamma
  xyz  = get.xyz(data, Gamma, px, py)
  rotx = rotmu[, 1:px, drop = FALSE]
  roty = rotmu[, px + 1:py, drop = FALSE]
  if (subcluster) {
    perms = unique.data.frame(gtools::permutations(n = G, r = G, v = rep(1:Gx, Gy), set = FALSE))
    perm.dists = apply(perms, 1, function(perm) {
      rotx |> split.data.frame(perm) |> sapply(\(x) sum(sqrt(colSums((t(x) - colMeans(x))^2)))) |> sum()
    })
    setx = as.integer(perms[which.min(perm.dists),])
    z_ngh = mc$z
  } else {
    mcx = Mclust(xyz$X, Gx, verbose = FALSE)$z
    mcy = Mclust(xyz$Y, Gy[1], verbose = FALSE)$z
    z_ngh = NULL
    for (g in 1:Gx) {z_ngh = cbind(z_ngh, zx[,g] * zy)}
  }

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster, diagonal = diagonal,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state, diagonal = diagonal)
  state$psi = apply(xyz$Z, 2, var)
  state = update.bars(state)

  for (i in 1:100) {
    state = e.step(state)
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }

  return(state)
}

init.random = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE) {
  if (px < 1) stop("px must be >= 1")
  if (py < 1) stop("py must be >= 1")

  n = nrow(data)
  p = ncol(data)
  pz = p - px - py

  Gamma = pracma::randortho(p)
  xyz  = get.xyz(data, Gamma, px, py)

  if (!subcluster) Gy = Gy[1]
  if (length(Gy) == 1) Gy = rep(Gy, Gx)
  G = sum(Gy)

  setx = rep(1:Gx, Gy)
  sety = if (subcluster) 1:G else rep(1:Gy[1], Gx)

  z_ngh = prop.table(matrix(runif(n * G), n, G), 1)

  tau.g.h = numeric(G)
  for (i in 1:G) {
    tau.g.h[i] =
      sum(z_ngh[,setx == setx[i]]) * sum(z_ngh[,sety == sety[i]])
  }
  tau.g.h = tau.g.h / sum(tau.g.h)

  state = list(
    data = data, px = px, py = py, Gx = Gx, Gy = Gy, G = G,
    regression = regression, subcluster = subcluster, diagonal = diagonal,
    setx = setx, sety = sety, Gamma = Gamma,
    X = xyz$X, Y = xyz$Y, Z = xyz$Z,
    mu.g = matrix(0, px, Gx),
    sigma.g = array(0, c(px, px, Gx)),
    beta.g.h = array(0, c(px + 1, py, max(sety))),
    sigma.g.h = array(0, c(py, py, max(sety))),
    psi = rep(0, pz),
    tau.g.h = tau.g.h,
    z_ngh = z_ngh
  )

  state = m.step(state, diagonal = diagonal)
  state$psi = state$psi * 0 + 10
  state = update.bars(state)

  for (i in 1:100) {
    state = e.step(state)
    state$Gamma = m_step_gamma_dual(
      state$G,
      state$sigma.bar.g.h.inv,
      state$mu.bar.g.h,
      state$z_ngh,
      state$data,
      state$Gamma
    )
  }
  return(state)
}

# Fitter Functions --------------------------------------------------------

fit.subclust = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE,
                        niter = 1000, mini.EM = 100, mini.EM.iter = 100, subsample = 0.5,
                        rng.seed = 12345, gamma.updates = 100, gamma.epsilon = 1e-8, init = list(init.mclust)) {
  state = NULL
  if (class(init) == "function") init = list(init)
  if (mini.EM > 0) {
    set.seed(rng.seed)
    old.seed = .Random.seed # capture rng state
    message("Mini-EM")
    mini.EM.par = expand.grid(rep = 1:mini.EM, init.fn = init)
    starts = pbapply::pbapply(mini.EM.par, 1, function(par) {
      state = tryCatch({
        sub.data = data[sort(sample(nrow(data), nrow(data) * subsample)), ]
        state = par[["init.fn"]](sub.data, px, py, Gx, Gy, regression, subcluster, diagonal)
        state = continue_em(state, mini.EM.iter, FALSE, FALSE, FALSE, diagonal, gamma.updates, gamma.epsilon)
        state$data = data
        state = update.xyz(state)
        state = e.step(state)
        state = continue_em(state, mini.EM.iter, FALSE, FALSE, FALSE, diagonal, gamma.updates, gamma.epsilon)
        ll = sum(logden.loglik(logden(state)))
        state$data = NULL # save some memory
        state$ll = ll
        state
      }, error = function(e) {
        print(e)
        list(ll = -Inf)
      })
      state
    }, simplify = FALSE)
    .Random.seed = old.seed # restore old rng state
    lls = sapply(starts, `[[`, "ll")

    state = starts[[which.max(lls)]]
    state$data = data
    message(sum(lls == -Inf), " failed during Mini-EM.")
    if (all(lls == -Inf)) stop("All Mini-EM failed.")
  } else {
    state = init[[1]](data, px, py, Gx, Gy, regression, subcluster, diagonal)
  }
  message("Actual Fit")
  state = continue_em(state, niter, FALSE, TRUE, TRUE, diagonal, gamma.updates, gamma.epsilon)
  return(state)
}

fit.subclust = function(data, px, py, Gx, Gy, regression = TRUE, subcluster = TRUE, diagonal = FALSE,
                        niter = 1000, mini.EM = 100, mini.EM.iter = 100, subsample = 0.5,
                        rng.seed = 12345, gamma.updates = 100, gamma.epsilon = 1e-8, init = list(init.mclust)) {
  state = NULL
  if (class(init) == "function") init = list(init)
  if (mini.EM > 0) {
    set.seed(rng.seed)
    old.seed = .Random.seed # capture rng state
    message("Mini-EM")
    mini.EM.par = expand.grid(rep = 1:mini.EM, init.fn = init)
    starts = pbapply::pbapply(mini.EM.par, 1, function(par) {
      state = tryCatch({
        sub.data = data[sort(sample(nrow(data), nrow(data) * subsample)), ]
        state = par[["init.fn"]](sub.data, px, py, Gx, Gy, regression, subcluster, diagonal)
        state = continue_em(state, mini.EM.iter, FALSE, FALSE, FALSE, diagonal, gamma.updates, gamma.epsilon)
        state$data = data
        state = update.xyz(state)
        state = e.step(state)
        state = continue_em(state, mini.EM.iter, FALSE, FALSE, FALSE, diagonal, gamma.updates, gamma.epsilon)
        ll = sum(logden.loglik(logden(state)))
        state$data = NULL # save some memory
        state$ll = ll
        state
      }, error = function(e) {
        print(e)
        list(ll = -Inf)
      })
      state
    }, simplify = FALSE)
    .Random.seed = old.seed # restore old rng state
    lls = sapply(starts, `[[`, "ll")

    state = starts[[which.max(lls)]]
    state$data = data
    message(sum(lls == -Inf), " failed during Mini-EM.")
    if (all(lls == -Inf)) stop("All Mini-EM failed.")
  } else {
    state = init[[1]](data, px, py, Gx, Gy, regression, subcluster, diagonal)
  }
  message("Actual Fit")
  state = continue_em(state, niter, FALSE, TRUE, TRUE, diagonal, gamma.updates, gamma.epsilon)
  return(state)
}

continue.em = function(state, niter, hold.z = FALSE, progress = TRUE, loglik = TRUE, diagonal = FALSE,
                       gamma.updates = 100, gamma.epsilon = 1e-8, loglik.eps = 1e-4) {
  if (progress) pb = txtProgressBar(0, niter, style = 3)
  ll = numeric(niter) * NA
  for(i in 1:niter) {
    if (!hold.z) state = e.step(state)
    state = m.step(state, diagonal = diagonal)
    state$Gamma = m_step_gamma_dual(
      state$G, state$sigma.bar.g.h.inv, state$mu.bar.g.h,
      state$z_ngh, state$data, state$Gamma,
      max_update = gamma.updates, epsilon = gamma.epsilon
    )
    state = update.xyz(state)
    if (progress) setTxtProgressBar(pb, i)
    if (loglik) {
      ll[i] = sum(logden.loglik(logden(state)))
      if (i > 1 && ll[i] - ll[i - 1] < loglik.eps) break
    }
  }
  if (loglik) state$ll = c(state$ll, ll)
  if (progress) close(pb)
  return(state)
}

continue.em.semi.sup = function(state, niter, hold.z = c(),
                                progress = TRUE, loglik = TRUE,
                                gamma.updates = 100, gamma.epsilon = 1e-8,
                                loglik.eps = 1e-4, diagonal = FALSE) {
  if (progress) pb = txtProgressBar(0, niter, style = 3)
  ll = numeric(niter) * NA
  z.held = state$z_ngh[hold.z, ]
  for(i in 1:niter) {
    state = e.step(state)
    state$z_ngh[hold.z, ] = z.held
    state = m.step(state, diagonal = diagonal)
    state$Gamma = m_step_gamma_dual(
      state$G, state$sigma.bar.g.h.inv, state$mu.bar.g.h,
      state$z_ngh, state$data, state$Gamma,
      max_update = gamma.updates, epsilon = gamma.epsilon
    )
    state = update.xyz(state)
    if (progress) setTxtProgressBar(pb, i)
    if (loglik) {
      ll[i] = sum(logden.loglik(logden(state)))
      if (i > 1 && ll[i] - ll[i - 1] < loglik.eps) break
    }
  }
  if (loglik) state$ll = c(state$ll, ll)
  if (progress) close(pb)
  return(state)
}

hold.z.ify = function(state, z.class) {
  for (i in 1:state$Gx) {state$z_ngh[z.class == i, state$setx != i] = 0}
  state$z_ngh = state$z_ngh / rowsums(state$z_ngh)
  return(state)
}

continue.em.semi.sup.z = function(state, niter, z.class, diagonal = FALSE,
                                  progress = TRUE, loglik = TRUE,
                                  gamma.updates = 100, gamma.epsilon = 1e-8,
                                  loglik.eps = 1e-4) {
  if (progress) pb = txtProgressBar(0, niter, style = 3)
  ll = numeric(niter) * NA
  for(i in 1:niter) {
    state = e.step(state)
    state = hold.z.ify(state, z.class)
    state = m.step(state, diagonal = diagonal)
    state$Gamma = m_step_gamma_dual(
      state$G, state$sigma.bar.g.h.inv, state$mu.bar.g.h,
      state$z_ngh, state$data, state$Gamma,
      max_update = gamma.updates, epsilon = gamma.epsilon
    )
    state = update.xyz(state)
    if (progress) setTxtProgressBar(pb, i)
    if (loglik) {
      ll[i] = sum(logden.loglik(logden(state)))
      if (i > 1 && ll[i] - ll[i - 1] < loglik.eps) break
    }
  }
  if (loglik) state$ll = c(state$ll, ll)
  if (progress) close(pb)
  return(state)
}

# Log-Likelihood ----------------------------------------------------------

logden = function(state) {
  # p = nrow(state$Gamma)
  # rot.data = cbind(state$X, state$Y, state$Z)
  #
  # zlog = vapply(1:state$G, function(k) {
  #   rot.data.center = t(t(rot.data) - state$mu.bar.g.h[[k]])
  #   siginv = state$sigma.bar.g.h.inv[[k]]
  #   mal = rowSums(rot.data.center %*% siginv * rot.data.center)
  #   -0.5 * mal + 0.5 * log(det(siginv)) + log(state$tau.g.h[k])
  # }, numeric(nrow(rot.data)))
  #
  # return(zlog - p/2*log(2*pi))

  zlog = logden_helper(state$G,
                       state$sigma.bar.g.h.inv, state$mu.bar.g.h,
                       state$tau.g.h,
                       state$X, state$Y, state$Z,
                       state$data)
  return(zlog)
}

mixlogden = function(state) {
  # pxy = state$px + state$py
  # rot.data = cbind(state$X, state$Y)
  #
  # zlog = vapply(1:state$G, function(k) {
  #   rot.data.center = t(t(rot.data) - state$mu.bar.g.h[[k]][1:pxy])
  #   siginv = state$sigma.bar.g.h.inv[[k]][1:pxy, 1:pxy]
  #   mal = rowSums(rot.data.center %*% siginv * rot.data.center)
  #   -0.5 * mal + 0.5 * log(det(siginv)) + log(state$tau.g.h[k])
  #   # constant -pxy/2 * log(2*pi) omitted
  # }, numeric(nrow(rot.data)))

  zlog = mixlogden_helper(state$px, state$py, state$G,
                          state$sigma.bar.g.h.inv, state$mu.bar.g.h,
                          state$tau.g.h,
                          state$X, state$Y)
  return(zlog)
}

logden.wts    <- function(z = NULL) {
  expz = exp(z)
  expz / rowSums(expz)
}
logden.wts = logden_wts

logden.loglik <- function(z = NULL) {
  log(rowSums(exp(z)))
}
logden.loglik = logden_loglik

# Expectation Step --------------------------------------------------------

e.step = function(state) {
  temp = logden.wts(mixlogden(state))
  # if (any(is.nan(temp))) return(state)
  state$z_ngh = temp
  return(state)
}

e.step2 = function(state) {
  state$z_ngh = e_step_helper(state$px, state$py, state$G,
                              state$sigma.bar.g.h.inv, state$mu.bar.g.h,
                              state$tau.g.h,
                              state$X, state$Y)
  return(state)
}

# Maximization Step -------------------------------------------------------

m.step = function(state, diagonal = FALSE) {
  out = m_step_primary(state$z_ngh, state$setx, state$Gx, state$X, diagonal = diagonal)
  pi.g = out$pi.g
  state$mu.g = out$mu.g
  state$sigma.g = out$sigma.g

  out = m_step_secondary(state$z_ngh, state$sety, state$Gy,
                         state$X, state$Y,
                         state$regression, state$subcluster, diagonal = diagonal)
  pi.g.h = out$pi.g.h
  state$sigma.g.h = out$sigma.g.h
  state$beta.g.h = out$beta.g.h

  state$psi = colMeans((state$Z)^2)

  pi.g = pi.g / nrow(state$z_ngh)
  pi.g.h = pi.g.h / nrow(state$z_ngh)
  for (i in 1:state$G) state$tau.g.h[i] = pi.g[state$setx[i]] * pi.g.h[state$sety[i]]

  state = update.bars(state)

  return(state)
}

m.step.gamma = function(state, updates = 15, epsilon = 1e-10) {
  state$Gamma = m_step_gamma(
    state$G,
    state$sigma.bar.g.h.inv,
    state$mu.bar.g.h,
    state$z_ngh,
    state$data,
    state$Gamma,
    max_update = updates,
    epsilon = epsilon
  )

  state = update.xyz(state)
  return(state)
}

m.step.gamma.tr = function(state, updates = 15, epsilon = 1e-10) {
  state$Gamma = m_step_gamma_tr(
    state$G,
    state$sigma.bar.g.h.inv,
    state$mu.bar.g.h,
    state$z_ngh,
    state$data,
    state$Gamma,
    max_update = updates,
    epsilon = epsilon
  )

  state = update.xyz(state)
  return(state)
}
