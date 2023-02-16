
gridSample <- function(xy, r, n=1, chess='') {
  
  if (inherits(xy, 'sf')) {
    xy <- st_coordinates(xy) #extracts coordinates from points
  }
  
  cell <- cellFromXY(r, xy) #extracts the cell id for each coordinate point
  uc <- unique(stats::na.omit(cell)) #extract NAs, points outside the grid, and keep only unique cells
  
  xy <- cbind(xy, cell, runif(nrow(xy)))#combine point coordinates, the cell id that contains them, and generate a random likelihood
  xy <- stats::na.omit(xy)
  xy <- unique(xy)
  
  
  xy <-  xy[order(xy[,4]), ] #order the points by their likelihood
  pts <- matrix(nrow=0, ncol=2)
  for (u in uc) { #for each unique cell with points
    ss <- subset(xy, xy[,3] == u) #get the points in that cell
    pts <- rbind(pts, ss[1:min(n, nrow(ss)), 1:2])#add the coordinates for the number of points allowed for each cell resolution
  }
  return(pts)
  
}

