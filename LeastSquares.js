// JS for website to find regression equation of a dataset
// Authors: Tyler Headley, Daniel Smith, Jeremy Tan

main();

// returns matrix AB, which is the matrix multiple of matrix A and B
function matrixMult(A, B) {

    if (A[0].length != B.length) {
        return "ERROR: matrices dimensions do not match";
    }

    // if A has dimensions a x b and B has dimensions b x c then new matrix is a x c
    const AB = new Array(A.length); 
    for (i = 0; i < AB.length; i++) {
        let row = new Array(B[0].length); 

        // calculate value at index i, j of the matrix by computing a dot product of A's ith row with B's jth column
        for (j = 0; j < row.length; j++) {
            let dotProduct = 0; 
            for (k = 0; k < B.length; k++) {
                dotProduct += A[i][k] * B[k][j];
            }
            row[j] = dotProduct;
        }
        AB[i] = row;
    }
    return AB;
}

// returns the transpose of matrix A, which is switching the rows with the columns
function transpose(A) {

    let AT = new Array(A[0].length); // flip dimensions of matrix
    for (let i = 0; i < AT.length; i++) {
        AT[i] = new Array(A.length); // create rows of transposed matrix
    }

    // iterate over new matrix switching indices from old matrix
    for (let i = 0; i < AT.length; i++) {
        for (let j = 0; j < AT[0].length; j++) {
            AT[i][j] = A[j][i];
        }
    }

    return AT;
}

// returns the determinant of the matrix A
function det(A) {

    // base case
    if (A.length == 2) {
        return A[0][0]*A[1][1] - A[0][1]*A[1][0]; // formula for det of a 2x2 matrix
    }

    else {
        let d = 0; 
        for (let i = 0; i < A.length; i++) {
            d += A[0][i] * cofactor(A, 0, i); // cofactor expansion along the top row
        }

        return d;
    }
}

// returns the cofactor of A at position x, y
function cofactor(A, x, y) {
    
    let c; // cofactor

    // base case
    if (A.length == 2) {
        c = A[-x+1][-y+1]; // element in opposite corner of 2x2 matrix
    }

    else {

        // shrink down the matrix to n-1 x n-1, excluding the xth row and yth column
        let a = new Array(A.length - 1); 
        for (let i = 0; i < a.length; i++) { 
            let arow = new Array(a.length);
            for (let j = 0; j < arow.length; j++) {

                
                let s = i, t = j;
                if (i >= x) {
                    s++; // skip over xth row
                }
                if (j >= y) {
                    t++; // skip over the yth column
                }

                arow[j] = A[s][t]; // add to new matrix
            }
            a[i] = arow;
        }
        c = det(a); // recursive call (cofactor and determinant call each other)
    }

    // if odd sum of indices, make cofactor negative
    if ((x+y) % 2 == 1) {
        c *= -1;
    }

    return c;
} 

// computes the inverse of the input matrix A
// known issue: because javascript store 10 decimal places, determinants near zero will be rounded down and produce an error
function inverse(A) {

    // inverse formula: A-1 = (1/|A|) * Adj(A)
    let inv = adj(A);
    let d = det(A);

    if (d == 0) {
        return "ERROR: matrix is not invertible since det = 0";
    }

    // multiply adjugate matrix by 1/det(A)
    for (let i = 0; i < A.length; i++) {
        for (let j = 0; j < A.length; j++) {
            inv[i][j] = inv[i][j] * (1/d);
        }
    }

    return inv;
}

// returns the adjugate of matrix A
// adjugate matrix is the transpose of its cofactor matrix 
// cofactor matrix is the matrix of cofactors of each element
function adj(A) {
    
    let cofactorMatrix = new Array(A.length); 

    // iterate over A
    for (let i = 0; i < A.length; i++) {
        let row = new Array(A[0].length);
        for (let j = 0; j < A[0].length; j++) {
            row[j] = cofactor(A,i,j); // find cofactor of element
        }
        cofactorMatrix[i] = row;
    }

    return transpose(cofactorMatrix);
}

// rounds each element of matrix A to the nearest 10000th decimal place
function roundMatrix(A) {

    for (let i = 0; i < A.length; i++) {
        for (let j = 0; j < A[i].length; j++) {
            A[i][j] = Math.round(A[i][j]*10000) / 10000; // Math.round() rounds to nearest integer
        }
    }
    return A;
}

// returns the coefficient vector for a polynomial least squares equation using the formula: x = ((A^TA)^-1)A^Tb
function leastSquares(data, degree) {

    // turn x coordinates into the matrix A for least squares
    // A has dimensions n x m, where n is the number of data points and m is the desired degree + 1
    const A = new Array(data.length);
    for (let i = 0; i < A.length; i++) {
        const row = new Array(degree+1);
        for (let j = 0; j < row.length; j++) {
            row[j] = Math.pow(data[i][0], j);
        }
        A[i] = row;
    }

    // create n x 1 column vector from y values
    b = new Array(data.length);
    for (let i = 0; i < b.length; i++) {
        b[i] = new Array(1);
        b[i][0] = data[i][1];
    }

    let AT = transpose(A);

    let ATA = matrixMult(AT, A);

    // known issue: if degree > number of data points there will be infinite solutions and this method won't work
    // known issue: if we have n data points and k have the same x value, if k > n - degree matrix ATA is not invertible 
        // potential fix: desmos seems to average the y values of the data points and consolidate into 1 data point
    let inv = inverse(ATA);

    let C = matrixMult(inv, AT);

    let x = matrixMult(C, b);

    return (roundMatrix(x)).map(arr => arr[0]); // turn n x 1 2d array into array of lenth n
}

// gets data from webpage, calls leastSquares and displays result 
function main() {

    const data = document.getElementById("data-input").value;
    const degree = Number(document.getElementById("degree").value);

    // regex creates a list of all matches of the form ([number],[number])
    // cannot have decimals as .xxxx; must be 0.xxxx; can be anything outside of parentheses and whitespace after comma
    let dataArray = data.match(/\(-?\d+\.?\d*,\s?-?\d+\.?\d*\)/g);
    

    for (let i = 0; i < dataArray.length; i++) {
        dataArray[i] = dataArray[i].match(/-?\d+\.?\d*/g); // extract numbers from order pair of coordinates

        // cast string to number
        dataArray[i][0] = Number(dataArray[i][0]);
        dataArray[i][1] = Number(dataArray[i][1]);
    }
    

    let coefficients = leastSquares(dataArray, degree);

    // use least squares coefficients to create polynomial
    let poly = "y = ";
    for (let i = 0; i < coefficients.length; i++) {
        if (i == 0) {
            poly += String(coefficients[i]);
        }
        if (i == 1) {
            poly += ` + ${coefficients[i]}x`;
        }
        if (i >= 2) {
            poly += ` + ${coefficients[i]}x<sup>${i}</sup>`;
        }
    }
        
    document.getElementById("equation").innerHTML = `Your regression equation is ${poly}`;
    
    graphData(dataArray, coefficients);
}

// displays both the input data and output polynomial regression
function graphData(dataArr, equation) {

    // Select the container element using d3 and normally
    const container = d3.select("#graph-container");
    const element = document.getElementById("graph-container");

    // remove any old data to refresh graph when new data is inputted
    while (element.firstChild) {
        element.removeChild(element.firstChild);
    }

    // Define the data
    const data = dataArr.map(point => ({ x: point[0], y: point[1] }));
    
    // find min and max x and y values using search algorithm
    let xMax = dataArr[0][0];
    let xMin = dataArr[0][0];
    let yMax = dataArr[0][1];
    let yMin = dataArr[0][1];

    for (let i = 1; i < dataArr.length; i++) {
        if (xMax < dataArr[i][0]){
            xMax = dataArr[i][0];
        }
        if (xMin > dataArr[i][0]){
            xMin = dataArr[i][0];
        }
        if (yMax < dataArr[i][1]){
            yMax = dataArr[i][1];
        }
        if (yMin > dataArr[i][1]){
            yMin = dataArr[i][1];
        }
    }

    const xRange = xMax - xMin;
    const yRange = yMax - yMin;

    // define the scales to display all data points, leaving some whitespace around the edges
    const xScale = d3.scaleLinear()
        .domain([xMin - .1*xRange, xMax + .1*xRange])
        .range([0, container.node().clientWidth]);

    const yScale = d3.scaleLinear()
        .domain([yMin - .1*yRange, yMax + .1*yRange])
        .range([container.node().clientHeight, 0]);

    // Create the scatterplot of input data
    container.append("svg")
        .attr("width", container.node().clientWidth)
        .attr("height", container.node().clientHeight)
        .selectAll("circle")
        .data(data)
        .enter()
        .append("circle")
        .attr("cx", d => xScale(d.x))
        .attr("cy", d => yScale(d.y))
        .attr("r", 4)
        .attr("fill", "blue");

    // graph polynomial

    // create scatterplot with enough points to create a solid line
    // known issue: for polynomials with very high slopes, points become separated
        // fix by calculating spacing of points dynamically based on slope
    const graph = d3.range(xMin - .1*xRange, xMax + .1*xRange, .0001*xRange);

    // convert 2d array to array of point objects
    for (let i = 0; i < graph.length; i++) {
        graph[i] = { x: graph[i], y: 0};
        for (let j = 0; j < equation.length; j++) {
            graph[i].y += equation[j] * graph[i].x ** j;
        }
    }

    // add to graph display
    container.append("svg")
        .attr("width", container.node().clientWidth)
        .attr("height", container.node().clientHeight)
        .selectAll("circle")
        .data(graph)
        .enter()
        .append("circle")
        .attr("cx", d => xScale(d.x))
        .attr("cy", d => yScale(d.y))
        .attr("r", 1)
        .attr("fill", "red");
    
}

