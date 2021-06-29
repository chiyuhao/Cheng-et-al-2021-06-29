draw.pairwise.venn.pro = function (area1, area2, cross.area, category = rep("", 2), euler.d = TRUE, 
                                   scaled = TRUE, inverted = FALSE, ext.text = TRUE, ext.percent = rep(0.05, 
                                                                                                       3), lwd = rep(2, 2), lty = rep("solid", 2), col = rep("black", 
                                                                                                                                                             2), fill = NULL, alpha = rep(0.5, 2), label.col = rep("black", 
                                                                                                                                                                                                                   3), cex = rep(1, 3), fontface = rep("plain", 3), fontfamily = rep("serif", 
                                                                                                                                                                                                                                                                                     3), cat.pos = c(-50, 50), cat.dist = rep(0.025, 2), 
                                   cat.cex = rep(1, 2), cat.col = rep("black", 2), cat.fontface = rep("plain", 
                                                                                                      2), cat.fontfamily = rep("serif", 2), cat.just = rep(list(c(0.5, 
                                                                                                                                                                  0.5)), 2), cat.default.pos = "outer", cat.prompts = FALSE, 
                                   ext.pos = rep(0, 2), ext.dist = rep(0, 2), ext.line.lty = "solid", 
                                   ext.length = rep(0.95, 2), ext.line.lwd = 1, rotation.degree = 0, 
                                   rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, 
                                   offset = 0, cex.prop = NULL, print.mode = "raw", sigdigs = 3, 
                                   set.maxcircle = F,
                                   maxcirclearea = 100,
                                   ...) 
{
  if (length(category) == 1) {
    category <- rep(category, 2)
  }
  else if (length(category) != 2) {
    flog.error("Unexpected parameter length for 'category'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'category'")
  }
  if (length(ext.percent) == 1) {
    ext.percent <- rep(ext.percent, 3)
  }
  else if (length(ext.percent) != 3) {
    flog.error("Unexpected parameter length for 'ext.percent'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'ext.percent'")
  }
  if (length(ext.pos) == 1) {
    ext.pos <- rep(ext.pos, 2)
  }
  else if (length(ext.pos) != 2) {
    flog.error("Unexpected parameter length for 'ext.pos'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'ext.pos'")
  }
  if (length(ext.dist) == 1) {
    ext.dist <- rep(ext.dist, 2)
  }
  else if (length(ext.dist) != 2) {
    flog.error("Unexpected parameter length for 'ext.dist'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'ext.dist'")
  }
  if (length(ext.length) == 1) {
    ext.length <- rep(ext.length, 2)
  }
  else if (length(ext.length) != 2) {
    flog.error("Unexpected parameter length for 'ext.length'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'ext.length'")
  }
  if (length(lwd) == 1) {
    lwd <- rep(lwd, 2)
  }
  else if (length(lwd) != 2) {
    flog.error("Unexpected parameter length for 'lwd'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'lwd'")
  }
  if (length(lty) == 1) {
    lty <- rep(lty, 2)
  }
  else if (length(lty) != 2) {
    flog.error("Unexpected parameter length for 'lty'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'lty'")
  }
  if (length(col) == 1) {
    col <- rep(col, 2)
  }
  else if (length(col) != 2) {
    flog.error("Unexpected parameter length for 'col'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'col'")
  }
  if (length(label.col) == 1) {
    label.col <- rep(label.col, 3)
  }
  else if (length(label.col) != 3) {
    flog.error("Unexpected parameter length for 'label.col'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'label.col'")
  }
  if (length(cex) == 1) {
    cex <- rep(cex, 3)
  }
  else if (length(cex) != 3) {
    flog.error("Unexpected parameter length for 'cex'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cex'")
  }
  if (length(fontface) == 1) {
    fontface <- rep(fontface, 3)
  }
  else if (length(fontface) != 3) {
    flog.error("Unexpected parameter length for 'fontface'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'fontface'")
  }
  if (length(fontfamily) == 1) {
    fontfamily <- rep(fontfamily, 3)
  }
  else if (length(fontfamily) != 3) {
    flog.error("Unexpected parameter length for 'fontfamily'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'fontfamily'")
  }
  if (length(fill) == 1) {
    fill <- rep(fill, 2)
  }
  else if (length(fill) != 2 & length(fill) != 0) {
    flog.error("Unexpected parameter length for 'fill'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'fill'")
  }
  if (length(alpha) == 1) {
    alpha <- rep(alpha, 2)
  }
  else if (length(alpha) != 2 & length(alpha) != 0) {
    flog.error("Unexpected parameter length for 'alpha'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'alpha'")
  }
  if (length(ext.line.lwd) != 1) {
    flog.error("Unexpected parameter length for 'ext.line.lwd'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'ext.line.lwd'")
  }
  if (length(cat.pos) == 1) {
    cat.pos <- rep(cat.pos, 2)
  }
  else if (length(cat.pos) != 2) {
    flog.error("Unexpected parameter length for 'cat.pos'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.pos'")
  }
  if (length(cat.dist) == 1) {
    cat.dist <- rep(cat.dist, 2)
  }
  else if (length(cat.dist) != 2) {
    flog.error("Unexpected parameter length for 'cat.dist'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.dist'")
  }
  if (length(cat.col) == 1) {
    cat.col <- rep(cat.col, 2)
  }
  else if (length(cat.col) != 2) {
    flog.error("Unexpected parameter length for 'cat.col'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.col'")
  }
  if (length(cat.cex) == 1) {
    cat.cex <- rep(cat.cex, 2)
  }
  else if (length(cat.cex) != 2) {
    flog.error("Unexpected parameter length for 'cat.cex'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.cex'")
  }
  if (length(cat.fontface) == 1) {
    cat.fontface <- rep(cat.fontface, 2)
  }
  else if (length(cat.fontface) != 2) {
    flog.error("Unexpected parameter length for 'cat.fontface'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.fontface'")
  }
  if (length(cat.fontfamily) == 1) {
    cat.fontfamily <- rep(cat.fontfamily, 2)
  }
  else if (length(cat.fontfamily) != 2) {
    flog.error("Unexpected parameter length for 'cat.fontfamily'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'cat.fontfamily'")
  }
  if (length(offset) != 1) {
    flog.error("Unexpected parameter length for 'offset'. Try using 'rotation.degree' to achieve non-vertical offsets", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter length for 'offset'. Try using 'rotation.degree' to achieve non-vertical offsets")
  }
  if (!(class(cat.just) == "list" & length(cat.just) == 2 & 
        length(cat.just[[1]]) == 2 & length(cat.just[[2]]) == 
        2)) {
    flog.error("Unexpected parameter format for 'cat.just'", 
               name = "VennDiagramLogger")
    stop("Unexpected parameter format for 'cat.just'")
  }
  if (!euler.d & scaled) {
    flog.error("Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.", 
               name = "VennDiagramLogger")
    stop("Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.")
  }
  if (offset > 1 | offset < 0) {
    flog.error("'Offset' must be between 0 and 1.  Try using 'rotation.degree = 180' to achieve offsets in the opposite direction.", 
               name = "VennDiagramLogger")
    stop("'Offset' must be between 0 and 1.  Try using 'rotation.degree = 180' to achieve offsets in the opposite direction.")
  }
  if (cross.area > area1 | cross.area > area2) {
    flog.error("Impossible: cross section area too large.", 
               name = "VennDiagramLogger")
    stop("Impossible: cross section area too large.")
  }
  cat.pos <- cat.pos + rotation.degree
  if (((cat.default.pos != "outer") & (cat.default.pos != 
                                       "text")) & cat.prompts) {
    flog.info("No default location recognized.  Automatically changing to 'outer'", 
              name = "VennDiagramLogger")
    cat.default.pos <- "outer"
  }
  if ((cat.default.pos == "outer") & cat.prompts) {
    flog.info("Placing category labels at default outer locations.  Use 'cat.pos' and 'cat.dist' to modify location.", 
              name = "VennDiagramLogger")
    flog.info(paste("Current 'cat.pos':", cat.pos[1], "degrees,", 
                    cat.pos[2], "degrees"), name = "VennDiagramLogger")
    flog.info(paste("Current 'cat.dist':", cat.dist[1], 
                    ",", cat.dist[2]), name = "VennDiagramLogger")
  }
  if ((cat.default.pos == "text") & cat.prompts) {
    flog.info("Placing category labels at default text locations.  Use 'cat.pos' and 'cat.dist' to modify location.", 
              name = "VennDiagramLogger")
    flog.info(paste("Current 'cat.pos':", cat.pos[1], "degrees,", 
                    cat.pos[2], "degrees"), name = "VennDiagramLogger")
    flog.info(paste("Current 'cat.dist':", cat.dist[1], 
                    ",", cat.dist[2]), name = "VennDiagramLogger")
  }
  max.circle.size = 0.2
  special.coincidental <- FALSE
  special.inclusion <- FALSE
  special.exclusion <- FALSE
  list.switch <- FALSE
  grob.list <- gList()
  if (!inverted) {
    tmp1 <- max(area1, area2)
    tmp2 <- min(area1, area2)
    if (tmp1 != area1) {
      list.switch <- TRUE
    }
    area1 <- tmp1
    area2 <- tmp2
    r1 <- sqrt(area1/pi)
    r2 <- sqrt(area2/pi)
    if (r2 == 0) {
      r2 <- 0.5 * r1
    }
    shrink.factor <- max.circle.size/r1
  }
  else {
    tmp1 <- max(area1, area2)
    tmp2 <- min(area1, area2)
    if (tmp1 != area1) {
      list.switch <- TRUE
    }
    area1 <- tmp1
    area2 <- tmp2
    r1 <- sqrt(area1/pi)
    r2 <- sqrt(area2/pi)
    if (r1 == 0) {
      r1 <- 0.5 * r2
    }
    shrink.factor <- max.circle.size/r2
  }
  
  if (set.maxcircle) {
    shrink.factor <- max.circle.size/sqrt(maxcirclearea/pi)
  }
  
  if (xor(list.switch, inverted)) {
    category <- rev(category)
    lwd <- rev(lwd)
    lty <- rev(lty)
    col <- rev(col)
    fill <- rev(fill)
    alpha <- rev(alpha)
    label.col <- rev(label.col)
    cex <- rev(cex)
    fontface <- rev(fontface)
    fontfamily <- rev(fontfamily)
    cat.pos <- rev(cat.pos)
    cat.dist <- rev(cat.dist)
    cat.col <- rev(cat.col)
    cat.cex <- rev(cat.cex)
    cat.fontface <- rev(cat.fontface)
    cat.fontfamily <- rev(cat.fontfamily)
    cat.just <- rev(cat.just)
    ext.pos <- rev(ext.pos)
    ext.length <- rev(ext.length)
  }
  r1 <- r1 * shrink.factor
  r2 <- r2 * shrink.factor
  if (area1 == area2 & area2 == cross.area) {
    special.coincidental <- TRUE
  }
  if (cross.area != 0 & (cross.area == area2 | cross.area == 
                         area1)) {
    special.inclusion <- TRUE
  }
  if (0 == cross.area) {
    special.exclusion <- TRUE
  }
  denom <- area1 + area2 - cross.area
  wrapLab <- function(num) {
    stri = ""
    if (print.mode[1] == "percent") {
      stri <- paste(signif(num * 100/denom, digits = sigdigs), 
                    "%", sep = "")
      if (isTRUE(print.mode[2] == "raw")) {
        stri <- paste(stri, "\n(", num, ")", sep = "")
      }
    }
    if (print.mode[1] == "raw") {
      stri <- num
      if (isTRUE(print.mode[2] == "percent")) {
        stri <- paste(stri, "\n(", paste(signif(num * 
                                                  100/denom, digits = sigdigs), "%)", sep = ""), 
                      sep = "")
      }
    }
    return(stri)
  }
  if (scaled & !special.inclusion & !special.exclusion & !special.coincidental) {
    d <- find.dist(area1, area2, cross.area, inverted = inverted)
    d <- d * shrink.factor
    x.centre.1 <- (1 + r1 - r2 - d)/2
    x.centre.2 <- x.centre.1 + d
    tmp <- VennDiagram::ellipse(x = x.centre.1, y = 0.5, 
                                a = ifelse(!inverted, r1, r2), b = ifelse(!inverted, 
                                                                          r1, r2), gp = gpar(lty = 0, fill = fill[1], 
                                                                                             alpha = alpha[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = x.centre.2, y = 0.5, 
                                a = ifelse(inverted, r1, r2), b = ifelse(inverted, 
                                                                         r1, r2), gp = gpar(lty = 0, fill = fill[2], 
                                                                                            alpha = alpha[2]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = x.centre.1, y = 0.5, 
                                a = ifelse(!inverted, r1, r2), b = ifelse(!inverted, 
                                                                          r1, r2), gp = gpar(lwd = lwd[1], lty = lty[1], 
                                                                                             col = col[1], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = x.centre.2, y = 0.5, 
                                a = ifelse(inverted, r1, r2), b = ifelse(inverted, 
                                                                         r1, r2), gp = gpar(lwd = lwd[2], lty = lty[2], 
                                                                                            col = col[2], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    if (length(cex.prop) > 0) {
      if (length(cex.prop) != 1) {
        flog.error("Value passed to cex.prop is not length 1", 
                   name = "VennDiagramLogger")
        stop("Value passed to cex.prop is not length 1")
      }
      func = cex.prop
      if (class(cex.prop) != "function") {
        if (cex.prop == "lin") {
          func = function(x) x
        }
        else if (cex.prop == "log10") {
          func = log10
        }
        else flog.error(paste0("Unknown value passed to cex.prop: ", 
                               cex.prop), name = "VennDiagramLogger")
        stop(paste0("Unknown value passed to cex.prop: ", 
                    cex.prop))
      }
      areas = c(area1 - cross.area, cross.area, area2 - 
                  cross.area)
      maxArea = max(areas)
      for (i in 1:length(areas)) {
        cex[i] = cex[i] * func(areas[i])/func(maxArea)
        if (cex[i] <= 0) 
          stop(paste0("Error in rescaling of area labels: the label of area ", 
                      i, " is less than or equal to zero"))
      }
    }
    if (ext.text) {
      area.1.pos <- x.centre.1 + ifelse(!inverted, -r1 + 
                                          ((2 * r1 - (r1 + r2 - d))/2), -r2 + ((2 * r2 - 
                                                                                  (r2 + r1 - d))/2))
      area.2.pos <- x.centre.2 + ifelse(!inverted, r2 - 
                                          ((2 * r2 - (r1 + r2 - d))/2), r1 - ((2 * r1 - 
                                                                                 (r2 + r1 - d))/2))
      if ((area1 - cross.area)/area1 > ext.percent[1] & 
          (area1 - cross.area)/area2 > ext.percent[1]) {
        tmp <- textGrob(label = wrapLab(ifelse(!inverted, 
                                               area1, area2) - cross.area), x = area.1.pos, 
                        y = 0.5, gp = gpar(col = label.col[1], cex = cex[1], 
                                           fontface = fontface[1], fontfamily = fontfamily[1]))
        grob.list <- gList(grob.list, tmp)
      }
      else {
        label.pos <- find.cat.pos(area.1.pos, 0.5, ext.pos[1], 
                                  ext.dist[1], r1)
        area.1.xpos <- label.pos$x
        area.1.ypos <- label.pos$y
        tmp <- textGrob(label = wrapLab(ifelse(!inverted, 
                                               area1, area2) - cross.area), x = area.1.xpos, 
                        y = area.1.ypos, gp = gpar(col = label.col[1], 
                                                   cex = cex[1], fontface = fontface[1], fontfamily = fontfamily[1]))
        grob.list <- gList(grob.list, tmp)
        tmp <- linesGrob(x = c(area.1.pos + ext.length[1] * 
                                 (area.1.xpos - area.1.pos), area.1.pos), y = c(0.5 + 
                                                                                  ext.length[1] * (area.1.ypos - 0.5), 0.5), 
                         gp = gpar(col = label.col[1], lwd = ext.line.lwd, 
                                   lty = ext.line.lty))
        grob.list <- gList(grob.list, tmp)
      }
      if ((area2 - cross.area)/area2 > ext.percent[2] & 
          (area2 - cross.area)/area1 > ext.percent[2]) {
        tmp <- textGrob(label = wrapLab(ifelse(inverted, 
                                               area1, area2) - cross.area), x = area.2.pos, 
                        y = 0.5, gp = gpar(col = label.col[3], cex = cex[3], 
                                           fontface = fontface[3], fontfamily = fontfamily[3]))
        grob.list <- gList(grob.list, tmp)
      }
      else {
        label.pos <- find.cat.pos(area.2.pos, 0.5, ext.pos[2], 
                                  ext.dist[2], r2)
        area.2.xpos <- label.pos$x
        area.2.ypos <- label.pos$y
        tmp <- textGrob(label = wrapLab(ifelse(inverted, 
                                               area1, area2) - cross.area), x = area.2.xpos, 
                        y = area.2.ypos, gp = gpar(col = label.col[3], 
                                                   cex = cex[3], fontface = fontface[3], fontfamily = fontfamily[3]))
        grob.list <- gList(grob.list, tmp)
        tmp <- linesGrob(x = c(area.2.pos + ext.length[1] * 
                                 (area.2.xpos - area.2.pos), area.2.pos), y = c(0.5 + 
                                                                                  ext.length[1] * (area.2.ypos - 0.5), 0.5), 
                         gp = gpar(col = label.col[3], lwd = ext.line.lwd, 
                                   lty = ext.line.lty))
        grob.list <- gList(grob.list, tmp)
      }
      if (cross.area/area2 > ext.percent[3] & cross.area/area1 > 
          ext.percent[3]) {
        tmp <- textGrob(label = wrapLab(cross.area), 
                        x = x.centre.1 + (d - ifelse(!inverted, r2, 
                                                     r1)) + (r1 + r2 - d)/2, y = 0.5, gp = gpar(col = label.col[2], 
                                                                                                cex = cex[2], fontface = fontface[2], fontfamily = fontfamily[2]))
        grob.list <- gList(grob.list, tmp)
      }
      else {
        cross.area.pos <- x.centre.1 + (d - r2) + (r1 + 
                                                     r2 - d)/2
        cross.pos <- find.cat.pos(cross.area.pos, 0.5, 
                                  ext.pos[1], ext.dist[1], r1 + r2)
        cross.area.xpos <- cross.pos$x
        cross.area.ypos <- cross.pos$y
        tmp <- textGrob(label = wrapLab(cross.area), 
                        x = cross.area.xpos, y = cross.area.ypos, 
                        gp = gpar(col = label.col[2], cex = cex[2], 
                                  fontface = fontface[2], fontfamily = fontfamily[2]))
        grob.list <- gList(grob.list, tmp)
        tmp <- linesGrob(x = c(cross.area.pos + ext.length[2] * 
                                 (cross.area.xpos - cross.area.pos), cross.area.pos), 
                         y = c(0.5 + ext.length[2] * (cross.area.ypos - 
                                                        0.5), 0.5), gp = gpar(col = label.col[2], 
                                                                              lwd = ext.line.lwd, lty = ext.line.lty))
        grob.list <- gList(grob.list, tmp)
      }
    }
    else {
      area.1.pos <- x.centre.1 + ifelse(!inverted, -r1 + 
                                          ((2 * r1 - (r1 + r2 - d))/2), -r2 + ((2 * r2 - 
                                                                                  (r2 + r1 - d))/2))
      tmp <- textGrob(label = wrapLab(ifelse(!inverted, 
                                             area1, area2) - cross.area), x = area.1.pos, 
                      y = 0.5, gp = gpar(col = label.col[1], cex = cex[1], 
                                         fontface = fontface[1], fontfamily = fontfamily[1]))
      grob.list <- gList(grob.list, tmp)
      area.2.pos <- x.centre.2 + ifelse(!inverted, r2 - 
                                          ((2 * r2 - (r1 + r2 - d))/2), r1 - ((2 * r1 - 
                                                                                 (r2 + r1 - d))/2))
      tmp <- textGrob(label = wrapLab(ifelse(inverted, 
                                             area1, area2) - cross.area), x = area.2.pos, 
                      y = 0.5, gp = gpar(col = label.col[3], cex = cex[3], 
                                         fontface = fontface[3], fontfamily = fontfamily[3]))
      grob.list <- gList(grob.list, tmp)
      tmp <- textGrob(label = wrapLab(cross.area), x = x.centre.1 + 
                        (d - ifelse(!inverted, r2, r1)) + (r1 + r2 - 
                                                             d)/2, y = 0.5, gp = gpar(col = label.col[2], 
                                                                                      cex = cex[2], fontface = fontface[2], fontfamily = fontfamily[2]))
      grob.list <- gList(grob.list, tmp)
    }
    if ("outer" == cat.default.pos) {
      cat.pos.1 <- find.cat.pos(x.centre.1, 0.5, (ifelse(!inverted, 
                                                         cat.pos[1], cat.pos[2]) + ifelse(xor(list.switch, 
                                                                                              inverted), 180, 0))%%360, cat.dist[1], ifelse(!inverted, 
                                                                                                                                            r1, r2))
      cat.pos.2 <- find.cat.pos(x.centre.2, 0.5, (ifelse(!inverted, 
                                                         cat.pos[2], cat.pos[1]) + ifelse(xor(list.switch, 
                                                                                              inverted), 180, 0))%%360, cat.dist[2], ifelse(!inverted, 
                                                                                                                                            r2, r1))
    }
    else if ("text" == cat.default.pos) {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], 
                                cat.dist[1])
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], 
                                cat.dist[2])
    }
    else {
      flog.error("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'", 
                 name = "VennDiagramLogger")
      stop("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'")
    }
    tmp <- textGrob(label = category[1], x = cat.pos.1$x, 
                    y = cat.pos.1$y, just = cat.just[[1]], gp = gpar(col = cat.col[1], 
                                                                     cex = cat.cex[1], fontface = cat.fontface[1], 
                                                                     fontfamily = cat.fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = category[2], x = cat.pos.2$x, 
                    y = cat.pos.2$y, just = cat.just[[2]], gp = gpar(col = cat.col[2], 
                                                                     cex = cat.cex[2], fontface = cat.fontface[2], 
                                                                     fontfamily = cat.fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
  }
  if (euler.d & special.inclusion & !special.coincidental) {
    if (inverted) {
      tmp1 <- area1
      tmp2 <- area2
      area1 <- tmp2
      area2 <- tmp1
    }
    if (!scaled & !inverted) {
      r1 <- 0.4
      r2 <- 0.2
    }
    if (!scaled & inverted) {
      r1 <- 0.2
      r2 <- 0.4
    }
    tmp <- VennDiagram::ellipse(x = 0.5, y = 0.5, a = r1, 
                                b = r1, gp = gpar(lty = 0, fill = fill[1], alpha = alpha[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.5 - offset * (r1 - 
                                                      r2), y = 0.5, a = r2, b = r2, gp = gpar(lty = 0, 
                                                                                              fill = fill[2], alpha = alpha[2]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.5, y = 0.5, a = r1, 
                                b = r1, gp = gpar(lwd = lwd[1], lty = lty[1], col = col[1], 
                                                  fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.5 - offset * (r1 - 
                                                      r2), y = 0.5, a = r2, b = r2, gp = gpar(lwd = lwd[2], 
                                                                                              lty = lty[2], col = col[2], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    area.2.pos <- 0.5 - offset * (r1 - r2)
    tmp <- textGrob(label = wrapLab(area2), x = area.2.pos, 
                    y = 0.5, gp = gpar(col = label.col[2], cex = cex[2], 
                                       fontface = fontface[2], fontfamily = fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
    if (!ext.text | !scaled) {
      area.1.pos <- (1 + r1 + r2 - offset * (r1 - r2))/2
      tmp <- textGrob(label = wrapLab(area1 - area2), 
                      x = area.1.pos, y = 0.5, gp = gpar(col = label.col[1], 
                                                         cex = cex[1], fontface = fontface[1], fontfamily = fontfamily[1]))
      grob.list <- gList(grob.list, tmp)
    }
    if (ext.text & scaled) {
      if (area2/area1 > 0.5) {
        area.1.pos <- (1 + r1 + r2 - offset * (r1 - 
                                                 r2))/2
        area.pos <- find.cat.pos(area.1.pos, 0.5, ext.pos[1], 
                                 ext.dist[1], r1)
        area.1.xpos <- area.pos$x
        area.1.ypos <- area.pos$y
        tmp <- textGrob(label = wrapLab(area1 - area2), 
                        x = area.1.xpos, y = area.1.ypos, gp = gpar(col = label.col[1], 
                                                                    cex = cex[1], fontface = fontface[1], fontfamily = fontfamily[1]))
        grob.list <- gList(grob.list, tmp)
        tmp <- linesGrob(x = c(area.1.pos + ext.length * 
                                 (area.1.xpos - area.1.pos), area.1.pos), y = c(0.5 + 
                                                                                  ext.length * (area.1.ypos - 0.5), 0.5), gp = gpar(col = label.col[1], 
                                                                                                                                    lwd = ext.line.lwd, lty = ext.line.lty))
        grob.list <- gList(grob.list, tmp)
      }
      else {
        area.1.pos <- (1 + r1 + r2 - offset * (r1 - 
                                                 r2))/2
        tmp <- textGrob(label = wrapLab(area1 - area2), 
                        x = area.1.pos, y = 0.5, gp = gpar(col = label.col[1], 
                                                           cex = cex[1], fontface = fontface[1], fontfamily = fontfamily[1]))
        grob.list <- gList(grob.list, tmp)
      }
    }
    if (cat.default.pos == "outer") {
      cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos[1], 
                                cat.dist[1], r1)
      cat.pos.2 <- find.cat.pos(0.5 - offset * (r1 - r2), 
                                0.5, cat.pos[2], cat.dist[2], r2)
    }
    else if (cat.default.pos == "text") {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], 
                                cat.dist[1])
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], 
                                cat.dist[2])
    }
    else {
      flog.error("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'", 
                 name = "VennDiagramLogger")
      stop("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'")
    }
    tmp <- textGrob(label = category[1], x = cat.pos.1$x, 
                    y = cat.pos.1$y, just = cat.just[[1]], gp = gpar(col = cat.col[1], 
                                                                     cex = cat.cex[1], fontface = cat.fontface[1], 
                                                                     fontfamily = cat.fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = category[2], x = cat.pos.2$x, 
                    y = cat.pos.2$y, just = cat.just[[2]], gp = gpar(col = cat.col[2], 
                                                                     cex = cat.cex[2], fontface = cat.fontface[2], 
                                                                     fontfamily = cat.fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
  }
  if (euler.d & special.coincidental) {
    tmp <- VennDiagram::ellipse(x = 0.5, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lty = 0, fill = fill[1], 
                                                               alpha = alpha[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.5, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lwd = lwd[1], lty = lty[1], 
                                                               col = col[1], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    area.1.pos <- 0.46
    tmp <- textGrob(label = wrapLab(area1), x = area.1.pos, 
                    y = 0.5, gp = gpar(col = label.col[2], cex = cex[2], 
                                       fontface = fontface[2], fontfamily = fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
    area.2.pos <- 0.54
    tmp <- textGrob(label = wrapLab(area2), x = area.2.pos, 
                    y = 0.5, gp = gpar(col = label.col[2], cex = cex[2], 
                                       fontface = fontface[2], fontfamily = fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = "(Coincidental)", x = 0.5, y = 0.45, 
                    gp = gpar(col = label.col[2], cex = cex[2], fontface = fontface[2], 
                              fontfamily = fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
    if (cat.default.pos == "outer") {
      cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos[1], 
                                cat.dist[1], max.circle.size)
      cat.pos.2 <- find.cat.pos(0.5, 0.5, cat.pos[2], 
                                cat.dist[2], max.circle.size)
    }
    else if (cat.default.pos == "text") {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], 
                                cat.dist[1])
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], 
                                cat.dist[2])
    }
    else {
      flog.error("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'", 
                 name = "VennDiagramLogger")
      stop("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'")
    }
    tmp <- textGrob(label = category[1], x = cat.pos.1$x, 
                    y = cat.pos.1$y, just = cat.just[[1]], gp = gpar(col = cat.col[1], 
                                                                     cex = cat.cex[1], fontface = cat.fontface[1], 
                                                                     fontfamily = cat.fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = category[2], x = cat.pos.2$x, 
                    y = cat.pos.2$y, just = cat.just[[2]], gp = gpar(col = cat.col[2], 
                                                                     cex = cat.cex[2], fontface = cat.fontface[2], 
                                                                     fontfamily = cat.fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
  }
  if (euler.d & special.exclusion) {
    if (!scaled) {
      r1 <- 0.2
      r2 <- 0.2
    }
    x.centre.1 <- (1 - 2 * (r1 + r2))/2 + r1 - sep.dist/2
    tmp <- VennDiagram::ellipse(x = x.centre.1, y = 0.5, 
                                a = r1, b = r1, gp = gpar(lty = 0, fill = fill[1], 
                                                          alpha = alpha[1]))
    grob.list <- gList(grob.list, tmp)
    x.centre.2 <- 1 - (1 - 2 * (r1 + r2))/2 - r2 + sep.dist/2
    tmp <- VennDiagram::ellipse(x = x.centre.2, y = 0.5, 
                                a = r2, b = r2, gp = gpar(lty = 0, fill = fill[2], 
                                                          alpha = alpha[2]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = x.centre.1, y = 0.5, 
                                a = r1, b = r1, gp = gpar(lwd = lwd[1], lty = lty[1], 
                                                          col = col[1], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = x.centre.2, y = 0.5, 
                                a = r2, b = r2, gp = gpar(lwd = lwd[2], lty = lty[2], 
                                                          col = col[2], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    area.1.pos <- x.centre.1
    tmp <- textGrob(label = wrapLab(area1), x = area.1.pos, 
                    y = 0.5, gp = gpar(col = label.col[1], cex = cex[1], 
                                       fontface = fontface[1], fontfamily = fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    area.2.pos <- x.centre.2
    tmp <- textGrob(label = wrapLab(area2), x = area.2.pos, 
                    y = 0.5, gp = gpar(col = label.col[3], cex = cex[3], 
                                       fontface = fontface[3], fontfamily = fontfamily[3]))
    grob.list <- gList(grob.list, tmp)
    if (cat.default.pos == "outer") {
      cat.pos.1 <- find.cat.pos(x.centre.1, 0.5, cat.pos[1], 
                                cat.dist[1], r1)
      cat.pos.2 <- find.cat.pos(x.centre.2, 0.5, cat.pos[2], 
                                cat.dist[2], r2)
    }
    else if (cat.default.pos == "text") {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], 
                                cat.dist[1])
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], 
                                cat.dist[2])
    }
    else {
      flog.error("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'", 
                 name = "VennDiagramLogger")
      stop("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'")
    }
    tmp <- textGrob(label = category[1], x = cat.pos.1$x, 
                    y = cat.pos.1$y, just = cat.just[[1]], gp = gpar(col = cat.col[1], 
                                                                     cex = cat.cex[1], fontface = cat.fontface[1], 
                                                                     fontfamily = cat.fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = category[2], x = cat.pos.2$x, 
                    y = cat.pos.2$y, just = cat.just[[2]], gp = gpar(col = cat.col[2], 
                                                                     cex = cat.cex[2], fontface = cat.fontface[2], 
                                                                     fontfamily = cat.fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
  }
  if ((!scaled & !euler.d) | (!scaled & euler.d & !special.inclusion & 
                              !special.exclusion & !special.coincidental)) {
    tmp <- VennDiagram::ellipse(x = 0.4, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lty = 0, fill = fill[1], 
                                                               alpha = alpha[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.6, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lty = 0, fill = fill[2], 
                                                               alpha = alpha[2]))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.4, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lwd = lwd[1], lty = lty[1], 
                                                               col = col[1], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    tmp <- VennDiagram::ellipse(x = 0.6, y = 0.5, a = max.circle.size, 
                                b = max.circle.size, gp = gpar(lwd = lwd[2], lty = lty[2], 
                                                               col = col[2], fill = "transparent"))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = wrapLab(area1 - cross.area), 
                    x = 0.3, y = 0.5, gp = gpar(col = label.col[1], 
                                                cex = cex[1], fontface = fontface[1], fontfamily = fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = wrapLab(area2 - cross.area), 
                    x = 0.7, y = 0.5, gp = gpar(col = label.col[3], 
                                                cex = cex[3], fontface = fontface[3], fontfamily = fontfamily[3]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = wrapLab(cross.area), x = 0.5, 
                    y = 0.5, gp = gpar(col = label.col[2], cex = cex[2], 
                                       fontface = fontface[2], fontfamily = fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
    if (cat.default.pos == "outer") {
      cat.pos.1 <- find.cat.pos(0.4, 0.5, cat.pos[1], 
                                cat.dist[1], max.circle.size)
      cat.pos.2 <- find.cat.pos(0.6, 0.5, cat.pos[2], 
                                cat.dist[2], max.circle.size)
    }
    else if (cat.default.pos == "text") {
      cat.pos.1 <- find.cat.pos(0.3, 0.5, cat.pos[1], 
                                cat.dist[1])
      cat.pos.2 <- find.cat.pos(0.7, 0.5, cat.pos[2], 
                                cat.dist[2])
    }
    else {
      flog.error("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'", 
                 name = "VennDiagramLogger")
      stop("Invalid value for 'cat.default.pos', should be either 'outer' or 'text'")
    }
    tmp <- textGrob(label = category[1], x = cat.pos.1$x, 
                    y = cat.pos.1$y, just = cat.just[[1]], gp = gpar(col = cat.col[1], 
                                                                     cex = cat.cex[1], fontface = cat.fontface[1], 
                                                                     fontfamily = cat.fontfamily[1]))
    grob.list <- gList(grob.list, tmp)
    tmp <- textGrob(label = category[2], x = cat.pos.2$x, 
                    y = cat.pos.2$y, just = cat.just[[2]], gp = gpar(col = cat.col[2], 
                                                                     cex = cat.cex[2], fontface = cat.fontface[2], 
                                                                     fontfamily = cat.fontfamily[2]))
    grob.list <- gList(grob.list, tmp)
  }
  
  #
  tmp <- VennDiagram::ellipse(x = 0.5, y = 0.5, a = max.circle.size,
                              b = max.circle.size, gp = gpar(lty = 0, fill = "transparent"))
  grob.list <- gList(grob.list, tmp)
  #
  
  grob.list <- adjust.venn(rotate.venn.degrees(grob.list, 
                                               rotation.degree, rotation.centre[1], rotation.centre[2]), 
                           ...)
  if (ind) {
    grid.draw(grob.list)
  }
  return(grob.list)
}
