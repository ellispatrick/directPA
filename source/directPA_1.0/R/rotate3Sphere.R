# spherical coordinates rotation

rotate3Sphere = function(T, direction = c(1,1,1)){
  
  # rotate to the direction that you want the test to face ie c(1,-1,-1)
  # T is a matrix of test statistics where the rows correspond to genes and the columns treatments
  
  a = direction
  b= c(1,1,1) # the original direction is c(1,1,1)
  
  if(sum(a == b)==3)return(T)
  if(sum(a == -b)==3)return(-T)
  
  a = a/sqrt(sum(a^2))
  b = b/sqrt(sum(b^2))
  u = c((a[2]*b[3] - a[3]*b[2]) , (a[3]*b[1] - a[1]*b[3]) ,(a[1]*b[2] - a[2]*b[1]) )
  u = u/sqrt(sum(u^2))
  R = matrix(0,3,3)
  the = acos(sum(a*b))
  R[1,1] = cos(the) + u[1]^2*(1-cos(the))
  R[2,2] = cos(the) + u[2]^2*(1-cos(the))
  R[3,3] = cos(the) + u[3]^2*(1-cos(the))           
  R[1,2] = u[1]*u[2]*(1-cos(the)) - u[3]*sin(the)                 
  R[1,3] = u[1]*u[3]*(1-cos(the)) + u[2]*sin(the)                 
  R[2,1] = u[2]*u[1]*(1-cos(the)) + u[3]*sin(the)                 
  R[2,3] = u[2]*u[3]*(1-cos(the)) - u[1]*sin(the)                 
  R[3,1] = u[3]*u[1]*(1-cos(the)) - u[2]*sin(the)                 
  R[3,2] = u[3]*u[2]*(1-cos(the)) + u[1]*sin(the)                 
  return(t(R%*%t(T)))
}
