#===============================================================================================
#
# grey2red
#
#===============================================================================================

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}


# Example of grey2red:

if (FALSE)
{
  par(mfrow = c(5,1))
  par(mar = c(1,3,1,1))
  n= 100
  barplot(rep(1, n), col = grey2red(n, 0, 1))
  barplot(rep(1, n), col = grey2red(n, 1, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 0.2))
  barplot(rep(1, n), col = grey2red(n, 0.5, 5.0))
}



#===============================================================================================
#
# Circle plot for generating VisANT-like plots
#
#===============================================================================================

circlePlot = function(
  adjacency,
  labels,
  order,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  lineColors = grey2red(50, 0.6, 1),
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  xMargin = 1-radii[1],
  yMargin = 1-radii[2],
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  
  ...)
{

  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);

  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = apply(abs(adjx), 2, sum)-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);

  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];

  actualCexPts = rep(0, n);
  for (node in 1:n)
  {
    cex = min.cex.points;
    if (variable.cex.points)
       cex = min.cex.points + (max.cex.points - min.cex.points) * 
                                  (oConn[node] - minConn)/(maxConn - minConn)
    actualCexPts[node] = cex
  }

  diag(oAdj) = 0;
  maxA = max(abs(oAdj));
  if (sum(oAdj < 0) > 0)
  {
     adjCol = numbers2colors(oAdj, signed = TRUE, lim = c(-maxA, maxA));
  } else {
     adjCol = numbers2colors(oAdj, signed = FALSE, lim = c(0, maxA));
  }


  ltA = oAdj;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;

  adjOrder = order(c(abs(ltA)))
  rows = row(oAdj)[adjOrder];
  cols = col(oAdj)[adjOrder];

  nLines = n*(n-1)/2;
  for (line in 1:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = oAdj[n1, n2];
    normA = abs(a)/maxA;

    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;

    #pRadius1 = par("cxy") * actualCexPts[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * actualCexPts[n2]/35;
    lineLen = sqrt( (x[n1] - x[n2])^2 + (y[n1] - y[n2])^2);
    x1 = x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen

    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }

  for (node in 1:n)
    points(x[node], y[node], pch = oPch[node], cex = actualCexPts[node], bg = oPBg[node], col = oPColors[node]);

  for (node in 1:n)
  {
    cex = min.cex.labels;
    if (variable.cex.labels)
       cex = min.cex.labels + (max.cex.labels - min.cex.labels) *
                                 (oConn[node] - minConn)/(maxConn - minConn)
    textWidth = strwidth(oLabs[node], cex = cex);
    textHeight = strheight(oLabs[node], cex = cex);
    if (variableLabelAngle)
    {
      ang = angles[node]/pi * 180;
      if (ang < 180) 
      {
         dir = 1;
      } else {
         dir = -1;
         ang = ang - 180;
      }
      ang = (90 - ang)/2
      xDir = 1;
      yDir = 1;
      cosAng = cos(ang/180*pi);
      sinAng = sin(ang/180*pi);
    } else {
      ang = 0;
      xDir = x[node];
      yDir = y[node];
      cosAng = 1;
      sinAng = 1;
      dir = 1;
    }
    angRad = ang/180*pi;
    pRadius = par("cxy") * actualCexPts[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(sum(c(cosAng^2, sinAng^2) * pRadius^2));
    rotMat = matrix( c(cosAng, sinAng, -sinAng, cosAng), 2, 2);
    labelShift = rotMat %*% as.matrix(c(textWidth, textHeight));
    text(x[node] + dir * xDir * (labelShift[1]/2 + cosAng * effPointRadius + xLabelOffset[node]), 
         y[node] + dir * yDir * (labelShift[2]/2 + sinAng * effPointRadius + yLabelOffset[node]), 
         labels = oLabs[node], adj = c(0.5, 0.5), 
         cex = cex, col = oLColors[node], srt = ang, xpd = TRUE);
  }

}

# Example of circle plot:

if (FALSE)
{

   sizeGrWindow(8,8)
   par(mfrow = c(1,1));
   nS = 100;
   nn = 30;
   mod = simulateModule(rnorm(nS), nn);
   adjacency = abs(cor(mod))^3;
   
   order = c(1:nn);
   
   labels = paste("Gene", c(1:nn));
   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, radii = c(0.6, 0.6));
   # Fix the overlaping labels - for now this requires semi-manual intervention.
   xOffset = rep(0.01, nn);
   xOffset[16] = -strwidth(labels[16])/3;
   circlePlot(adjacency, labels, order, xLabelOffset = xOffset);


   # Plot two circles in one plot

   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, center = c(-0.5, -0.5), 
              radii = c(0.35, 0.35));

   circlePlot(adjacency, labels, order, startNewPlot = FALSE, 
              variable.cex.labels = FALSE, center = c(0.5, 0.5), 
              radii = c(0.35, 0.35));

   
}
