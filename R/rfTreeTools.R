# Walk the tree depth-first and find the height of the tree
get.tree.height <- function(tree) {
  current.height <- 0
  max.height <- 0
  prior.index.stack <- c(0)
  current.index <- 1

  while (current.index != 0) {
    max.height <- max(max.height, current.height)

    if (all(tree[current.index, c(1,2)] == 0)) {
      #backtracking - pop the parent off the stack
      current.index <- prior.index.stack[1]
      prior.index.stack <- prior.index.stack[-1]

      if (current.index != 0) {
        # delete node we just visited
        if (tree[current.index, 1] != 0) {
          tree[current.index, 1] <- 0
        } else if (tree[current.index, 2] != 0) {
          tree[current.index, 2] <- 0
        }
      }

      #decrement our height
      current.height <- current.height - 1

    } else {
      #push this level onto the stack
      prior.index.stack <- c(current.index, prior.index.stack)
      #increment our height
      current.height <- current.height + 1
      if (tree[current.index, 1] != 0) {
        current.index <- tree[current.index, 1]
      } else if (tree[current.index, 2] != 0) {
        current.index <- tree[current.index, 2]
      }
    }
  }
  return(max.height)
}

navigate.rf.tree <- function(tree, data, height=NULL) {
  current.height <- 0
  current.index <- 1
  data <- as.data.frame(data)

  left.index <- 1
  right.index <- 2
  split.var.index <- 3
  split.point.index <- 4


  split.var <- tree[current.index, split.var.index]

  if (is.null(height)) {
    height <- nrow(tree)
  }

  while (current.height < height & split.var != 0) {
    split.var <- tree[current.index, split.var.index]
    if (split.var == 0) {
      # We've hit a leaf node - break out
      break
    }

    if (is.data.frame(data)) {
      data.split.var <- data[1, split.var]
    } else {
      data.split.var <- data[split.var]
    }

    if (is.numeric(data.split.var)) {
      # Continuous
      split.point <- tree[current.index, split.point.index]

      if (data.split.var <= split.point) {
        current.index <- tree[current.index, left.index]
      } else {
        current.index <- tree[current.index, right.index]
      }

    } else if (is.factor(data.split.var)) {
      # Assume categorical
      split.point <- as.numeric(intToBits(tree[current.index, split.point.index]))
      factor.level <- as.numeric(data.split.var)
      #print(split.point)
      #print(factor.level)
      if (split.point[factor.level] == 1) {
        current.index <- tree[current.index, left.index]
        #print("l")
      } else {
        current.index <- tree[current.index, right.index]
        #print("r")
      }

    } else {
      stop(paste("Variable", split.var, "is neither numeric nor a factor."))
    }

    current.height <- current.height + 1
  }

  return(current.index)
}

#' Create Additional Params List for Random Forest Kernel
#'
#' Takes a \code{\link[randomForest]{randomForest}} object and extracts the
#' components required by \code{\link{randomForestKernel}} and
#' \code{\link{randomForestKernelGrad}}
#'
#' @param rf a randomForest object
#'
#' @return a list containing the additional params for \code{\link{randomForestKernel}}.
#' @export
#'
#' @importFrom randomForest getTree
create.rf.additional.params <- function(rf) {
  additional.params <- list()
  forest <- list()
  num.trees <- rf$ntree
  heights <- numeric(num.trees)
  for (i in 1:num.trees) {
    forest[[i]] <- getTree(rf, i)
    heights[i] <- sample(0:getTreeHeight(forest[[i]]), 1)
  }
  is.factor.vec <- vector(mode="logical", length=length(rf$forest$xlevels))
  for (i in seq_along(rf$forest$xlevels)) {
    xlevels.it <- rf$forest$xlevels[[i]]
    if (is.character(xlevels.it)) {
      is.factor.vec[i] <- TRUE
    } else {
      is.factor.vec[i] <- FALSE
    }
  }
  additional.params$forest <- forest
  additional.params$heights <- heights
  additional.params$isFactor <- is.factor.vec

  return(additional.params)
}
