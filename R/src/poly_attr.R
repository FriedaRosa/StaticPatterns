
## ======================================= ##
#                                           #
#     Function 7: Spatial Geometries        #
#     Function Credits: XXX ANONYMIZED      #
#                                           #
## ======================================= ##


# Function to get the max elongation per polygon, north-south length or east-west length, circularity, etc...
poly_attr <- function(x, type = NULL) {
  loadNamespace("geosphere")
  loadNamespace("terra")
  loadNamespace("tidyterra")

  # Convert "sf" to SpatVector
  if (inherits(x, "sf") == T) {
    x <- terra::vect(x)
  }
  # Reproject to WGS84. Currently, "terra" calculate distance in meters when given latitude and longitude points.
  x <- project(x, "epsg:4326")
  # Get East-West distance (longitudinal extension)
  if (type == "ewDist") {
    vert <- crds(as.points(ext(x)), df = T)
    dimNames <- list(c("point"),c("x","y"))
    sw <- matrix(c(min(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    se <- matrix(c(max(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    res <- distance(sw, se, lonlat = T)[[1]]
  }
  # Get South-North distance (latitudinal extension)
  if (type == "nsDist") {
    vert <- crds(as.points(ext(x)), df = T)
    dimNames <- list(c("point"),c("x","y"))
    sw <- matrix(c(min(vert[["x"]]), min(vert[["y"]])), ncol = 2, dimnames = dimNames)
    nw <- matrix(c(min(vert[["x"]]), max(vert[["y"]])), ncol = 2, dimnames = dimNames)
    res <- distance(sw, nw, lonlat = T)[[1]]
  }
  # Get max distance between opposite vertices
  if (type == "maxDist") {
    vert <- crds(as.points(convHull(x)))
    res <- distance(vert, lonlat = T) %>% max()
  }
  # Get elongation ratio along longest axis
  if (type == "elonRatio") {
    convexHull <- convHull(x)
    vert <- crds(as.points(convexHull), df = T)
    dist <- as.data.frame.table(as.matrix(distance(vert, lonlat = T)), responseName = "distance") %>%
      slice_max(., distance)
    axisPoints <- vert[c(dist[[1]][[1]],dist[[2]][[1]]),]
    axisPoints <- arrange(axisPoints, desc(y))
    rotation <- -1*bearing(axisPoints[2,],axisPoints[1,])
    rotHull <- spin(convexHull, rotation)
    ext <- ext(convexHull)
    df <- as.vector(distance(crds(as.points(ext)), lonlat = T))
    df <- sort(df)
    length <- mean(df[[3]], df[[4]])
    width <- mean(df[[1]], df[[2]])
    res <- 1 - (width / length)
    # res <- list()
    # res[["vert"]] <- vert
    # res[["dist"]] <- dist
    # res[["axispoints"]] <- axisPoints
    # res[["bearing"]] <- rotation
    # res[["rotHull"]] <- rotHull
  }
  # Circularity
  if (type == "circ") {
    perimeter <- perim(x)
    area <- expanse(x)
    res <- (perimeter^2) / area
  }
  # Normalized circularity
  if (type == "circNorm") {
    perimeter <- perim(x)
    area <- expanse(x)
    res <- (perimeter^2) / (4 * pi * area)
  }
  # Major length of the minimum rectangle
  if (type == "lengthMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle), df = T), lonlat = T))
    df <- sort(df)
    res <- mean(df[[3]], df[[4]])
  }
  # Width of the minimum rectangle
  if (type == "widthMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle), df = T), lonlat = T))
    df <- sort(df)
    res <- mean(df[[1]], df[[2]])
  }
  # Elongation ratio of minimal encasing rectangle
  # (from here: Dražić, S., Ralević, N., & Žunić, J. (2010).
  # Shape elongation from optimal encasing rectangles.
  # Computers & Mathematics with Applications, 60(7), 2035–2042.
  # https://doi.org/10.1016/j.camwa.2010.07.043)
  if (type == "elonMinRect") {
    minRectangle <- minRect(x)
    df <- as.vector(distance(crds(as.points(minRectangle)), lonlat = T))
    df <- sort(df)
    length <- mean(df[[3]], df[[4]])
    width <- mean(df[[1]], df[[2]])
    res <- 1 - (width / length)
  }
  # Related circumscribing circle
  if (type == "relCirc") {
    # if (values(x)$datasetID[1] == 26) {
    #   x <- project(x, vars$crs[4])
    # }
    circle <- minCircle(x)
    areaCircle <- expanse(circle, transform=FALSE)
    area <- expanse(x)
    res <- 1-(area/areaCircle)
  }
  # Linearity index
  if (type == "lin") {
    hull <- convHull(x)
    df <- crds(as.points(hull), df = T)
    lm <- lm(y ~ x, data = df)
    res <- summary(lm)$r.squared
  }
  # North bearing of the minimum rectangle
  if (type == "bearingMinRect") {
    minRectangle <- minRect(x)
    df <- crds(as.points(minRectangle), df = T)
    cor <- cor(df[["x"]],df[["y"]])
    point1 <- slice_min(df, y)
    if (cor > 0){
      point2 <-slice_max(df, x)
    } else{
      point2 <-slice_min(df, x)
    }
    res <- bearing(point1, point2)
  }
  # Get bearing along longest axis
  if (type == "bearing") {
    convexHull <- convHull(x)
    vert <- crds(as.points(convexHull), df = T)
    dist <- as.data.frame.table(as.matrix(distance(vert, lonlat = T)), responseName = "distance") %>%
      slice_max(., distance)
    axisPoints <- vert[c(dist[[1]][[1]],dist[[2]][[1]]),]
    axisPoints <- arrange(axisPoints, desc(y))
    res <- bearing(axisPoints[2,],axisPoints[1,])
  }
  res <- as.numeric(res)
  return(res)
}



