# polar coordinates rotation

rotate2Sphere = function(T, degree = 0){
  # rotate to the direction that you want the test to face ie pi/4 degree
  # T is a matrix of test statistics where the rows correspond to genes and the columns treatments

  theta = atan2(T[,2], T[,1]) + degree  
  r = sqrt(T[,1]^2 + T[,2]^2)
  x2 = r*cos(theta)
  y2 = r*sin(theta)
  return(cbind(x2, y2))
}
