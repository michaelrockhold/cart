//
//  Cart.swift
//  swiftcart
//
//  Created by Michael Rockhold on 8/18/20.
//  Copyright © 2020 Michael Rockhold. All rights reserved.
//

import Foundation
import Accelerate

protocol Population {
    associatedtype ElementType
    
    var height: Int { get }
    var width: Int { get }
    var data: [ElementType] { get }
}

func clip(f: Float, min: Int, max: Int) -> Int {
    let i = Int(f)
    if i < min {
        return min
    }
    else if i > max {
        return max
    }
    else {
        return i
    }
}

func clipf(f: Float, min: Int, max: Int) -> Float {
    if f < Float(min) {
        return Float(min)
    }
    else if f > Float(max) {
        return Float(max)
    }
    else {
        return f
    }
}

struct Point<T> {
    var x: T
    var y: T
}

struct Matrix2D<T> {
    let rows: Int, columns: Int
    var data: [T]
    init(rows: Int, columns: Int, initialValue: T) {
        self.rows = rows
        self.columns = columns
        data = [T](repeating: initialValue, count: rows * columns)
    }
    
    init<P: Population>(fromPopulation pop: P) where P.ElementType == Float, T == Float {
        self.rows = pop.height
        self.columns = pop.width
        self.data = pop.data
    }
    
    func indexIsValid(row: Int, column: Int) -> Bool {
        return row >= 0 && row < rows && column >= 0 && column < columns
    }
    subscript(row: Int, column: Int) -> T {
        get {
            assert(indexIsValid(row: row, column: column), "Index out of range")
            return data[(row * columns) + column]
        }
        set {
            assert(indexIsValid(row: row, column: column), "Index out of range")
            data[(row * columns) + column] = newValue
        }
    }
    
    func rowRange(_ row: Int) -> Range<Int> {
        let start = columns * row
        let next = columns * (row + 1)
        return start..<next
    }
    
    typealias srcAction = (UnsafePointer<T>) -> Void
    typealias dstAction = (UnsafeMutablePointer<T>) -> Void
    
    func withUnsafePointerToRow(_ row: Int, action: srcAction) {
        self.data[rowRange(row)].withUnsafeBufferPointer { buffer in
            action( buffer.baseAddress! )
        }
    }
    
    mutating func withUnsafeMutablePointerToRow(_ row: Int, action: dstAction) {
        self.data[rowRange(row)].withUnsafeMutableBufferPointer { buffer in
            action( buffer.baseAddress! )
        }
    }
    
}

struct DCT2DPlan {
    enum Direction {
        case forward
        case backward
    }
    
    let direction: Direction
    let height: Int
    let width: Int
    let rowwiseSetup: vDSP.DCT
    let columnwiseSetup: vDSP.DCT
    
    let directionmap:[Direction:vDSP.DCTTransformType] = [
        .forward: vDSP.DCTTransformType.II,
        .backward: vDSP.DCTTransformType.III
    ]
    
    init(height: Int, width: Int, direction: Direction) {
        self.height = height
        self.width = width
        self.direction = direction
        
        self.rowwiseSetup = vDSP.DCT(count: width, transformType: directionmap[direction]!)!
        self.columnwiseSetup = vDSP.DCT(count: height, transformType: directionmap[direction]!)!
    }
    
    func execute(src: Matrix2D<Float>, dst: inout Matrix2D<Float>) {
        
//        for row in 0..<height {
//            src.withUnsafePointerToRow(row) { srcbuffer in
//                dst.withUnsafeMutablePointerToRow(row) { dstbuffer in
//                    rowwiseSetup.transform(srcbuffer, result: &dstbuffer)
//                }
//            }
//        }
//        
//        vDSP_mtrans(&dst.data, 1, &dst.data, 1, vDSP_Length(dst.columns), vDSP_Length(dst.rows))
//        
//        // after transposition, there are 'width' rows
//        for row in 0..<width {
//            dst.withUnsafePointerToRow(row) { srcbuffer in
//                dst.withUnsafeMutablePointerToRow(row) { dstbuffer in
//                    columnwiseSetup.transform(srcbuffer, result: &dstbuffer)
//                }
//            }
//        }
//        
//        vDSP_mtrans(&dst.data, 1, &dst.data, 1, vDSP_Length(dst.rows), vDSP_Length(dst.columns))
    }
}

struct Cart {
    
    let INITH: Float = 0.001          // Initial size of a time-step
    let TARGETERROR: Float = 0.01     // Desired accuracy per step in pixels
    let MAXRATIO: Float = 4.0         // Max ratio to increase step size by
    let EXPECTEDTIME: Float = 1.0e8   // Guess as to the time it will take, used to estimate completion
    let PI: Float = 3.1415926535897932384626
    
    let SNAPSHOTCOUNT = 5

    let WIDTH: Int
    let HEIGHT: Int

    let rhotplan: [DCT2DPlan]   // Plan for rho(t) back-transform at time t
    

    var rhot: [Matrix2D<Float>]   // Pop density at time t (SNAPSHOTCOUNT snaps needed)
    var fftrho: Matrix2D<Float>   // FT of initial density
    var fftexpt: Matrix2D<Float>  // FT of density at time t
    
    var vxt: [Matrix2D<Float>]  // x-velocity at time t
    var vyt: [Matrix2D<Float>]  // y-velocity at time t
    
    var expky: [Float]          // Array needed for the Gaussian convolution
    
    
    init<P: Population>(population: P) where P.ElementType == Float  {
        WIDTH = population.width
        HEIGHT = population.height
        let popData = Matrix2D<Float>(fromPopulation: population)
        
        rhot = [Matrix2D](repeating: Matrix2D(rows: HEIGHT, columns: WIDTH, initialValue: 0.0), count: SNAPSHOTCOUNT)
        
        fftrho = Matrix2D<Float>(rows: HEIGHT, columns: WIDTH, initialValue: 0.0)
        fftexpt = Matrix2D<Float>(rows: HEIGHT, columns: WIDTH, initialValue: 0.0)
        
        vxt = [Matrix2D<Float>](repeating: Matrix2D<Float>(rows: HEIGHT+1, columns: WIDTH+1, initialValue: 0.0), count: SNAPSHOTCOUNT)
        vyt = [Matrix2D<Float>](repeating: Matrix2D<Float>(rows: HEIGHT+1, columns: WIDTH+1, initialValue: 0.0), count: SNAPSHOTCOUNT)
        
        expky = [Float](repeating: 0.0, count: HEIGHT)
        
        /* Make plans for the back transforms;
         * You'll execute these later with src=fftexpt, dst=rhot[i]
         */
        
        rhotplan = [DCT2DPlan](repeating: DCT2DPlan(height: HEIGHT, width: WIDTH, direction: .backward), count:SNAPSHOTCOUNT)
        
        /* Function to calculate the discrete cosine transform of the input data.
         * assumes its input is an fftw_malloced array in column-major form with
         * size WIDTH*HEIGHT */
        
        DCT2DPlan(height: HEIGHT, width: WIDTH, direction: .forward).execute(src: popData, dst: &fftrho)
    }

    /* Function to calculate the population density at arbitrary time by back-
     * transforming and put the result in a particular rhot[] snapshot array.
     * Calculates unnormalized densities, since FFTW gives unnormalized back-
     * transforms, but this doesn't matter because the cartogram method is
     * insensitive to variation in the density by a multiplicative constant */
    
    mutating func density(t: Float, s: Int) {
        
        /* Calculate the expky array, to save time in the next part */
        
        for iy in 0..<HEIGHT {
            let ky = PI * Float(iy) / Float(HEIGHT)
            expky[iy] = exp(-ky*ky*t)
        }
        
        /* Multiply the FT of the density by the appropriate factors */
        
        for ix in 0..<WIDTH {
            let kx = PI * Float(ix) / Float(WIDTH)
            let expkx = exp(-kx*kx*t)
            for iy in 0..<HEIGHT {
                fftexpt[iy, ix] = expkx * expky[iy] * fftrho[iy, ix]
            }
        }
        
        /* Perform the back-transform */
        
        rhotplan[s].execute(src: fftexpt, dst: &rhot[s])
    }
    
    
    /* Function to calculate the velocity at all integer grid points for a
     * specified snapshot */
    
    mutating func vgrid(_ s: Int) {
        
        /* Do the corners */
        
        vxt[s][0, 0] = 0.0
        vyt[s][0, 0] = 0.0
        vxt[s][0, WIDTH] = 0.0
        vyt[s][0, WIDTH] = 0.0
        
        vxt[s][HEIGHT, 0] = 0.0
        vyt[s][HEIGHT, 0] = 0.0
        vxt[s][HEIGHT, WIDTH] = 0.0
        vyt[s][HEIGHT, WIDTH] = 0.0
        
        /* Do the top border */
        
        var r11: Float = rhot[s][0, 0]
        for c in 1..<WIDTH {
            let r01 = r11
            r11 = rhot[s][0, c]
            vxt[s][0, c] = -2*(r11-r01)/(r11+r01)
            vyt[s][0, c] = 0.0
        }
        
        /* Do the bottom border */
        
        var r10 = rhot[s][HEIGHT-1, 0]
        for c in 1..<WIDTH {
            let r00 = r10
            r10 = rhot[s][HEIGHT-1, c]
            vxt[s][HEIGHT, c] = -2*(r10-r00)/(r10+r00)
            vyt[s][HEIGHT, c] = 0.0
        }
        
        /* Left edge */
        
        r11 = rhot[s][0, 0]
        for r in 1..<HEIGHT {
            let r10 = r11
            r11 = rhot[s][r, 0]
            vxt[s][r, 0] = 0.0
            vyt[s][r, 0] = -2*(r11-r10)/(r11+r10)
        }
        
        /* Right edge */
        
        var r01 = rhot[s][0, WIDTH-1]
        for r in 1..<HEIGHT {
            let r00 = r01
            r01 = rhot[s][r, WIDTH-1]
            vxt[s][r, WIDTH] = 0.0
            vyt[s][r, WIDTH] = -2*(r01-r00)/(r01+r00)
        }
        
        /* Now do all the points in the middle */
        
        for c in 1..<WIDTH {
            var r01 = rhot[s][0, c-1]
            var r11 = rhot[s][0, c]
            for r in 1..<HEIGHT {
                let r00 = r01
                let r10 = r11
                r01 = rhot[s][r, c-1]
                r11 = rhot[s][r, c]
                let mid = r10 + r00 + r11 + r01
                vxt[s][r, c] = -2*(r10-r00+r11-r01)/mid
                vyt[s][r, c] = -2*(r01-r00+r11-r10)/mid
            }
        }
    }
    
    
    /* Function to calculate the velocity at an arbitrary point from the grid
     * velocities for a specified snapshot by interpolating between grid
     * points.  If the requested point is outside the boundaries, we
     * extrapolate (ensures smooth flow back in if we get outside by mistake,
     * although we should never actually do this because function cart_twosteps()
     * contains code to prevent it) */
    
    func velocity(_ rx: Float, _ ry: Float, _ s: Int) -> (Float, Float) {
                
        /* Deal with the boundary conditions */
        
        let col = clip(f: rx, min: 0, max: WIDTH-1)
        let row = clip(f: ry, min: 0, max: HEIGHT-1)
        
        /* Calculate the weights for the bilinear interpolation */
        
        let dx = rx - Float(col)
        let dy = ry - Float(row)
        
        let dx1m = 1.0 - dx
        let dy1m = 1.0 - dy
        
        let w11 = dx1m * dy1m
        let w21 = dx * dy1m
        let w12 = dx1m * dy
        let w22 = dx * dy
        
        /* Perform the interpolation for x and y components of velocity */
        let vxp = w11*vxt[s][row, col]
                + w21*vxt[s][row, col+1]
                + w12*vxt[s][row+1, col]
                + w22*vxt[s][row+1, col+1]
        let vyp = w11*vyt[s][row, col]
                + w21*vyt[s][row, col+1]
                + w12*vyt[s][row+1, col]
                + w22*vyt[s][row+1, col+1]
        return (vxp, vyp)
    }
    
    
    /* Function to integrate 2h time into the future two different ways using
     * four-order Runge-Kutta and compare the differences for the purposes of
     * the adaptive step size.  Parameters are:
     *   *pointx = array of x-coords of points
     *   *pointy = array of y-coords of points
     *   npoints = number of points
     *   t = current time, i.e., start time of these two steps
     *   h = delta t
     *   s = snapshot index of the initial time
     *   WIDTH, HEIGHT = size of grid
     * Returns (errorp, drp, spp) as a tuple:
     *   errorp = the maximum integration error found for any polygon vertex for
     *             the complete two-step process
     *   drp = maximum distance moved by any point
     *   spp = the snapshot index for the final function evaluation
     */
    
    mutating func twosteps(points: [Point<Float>],
                           t: Float, h: Float, s: Int) -> (Float, Float, Int, [Point<Float>]) {
        
        /* Do the big combined (2h) RK step */
        func bigRKstep(_ rx1: Float, _ ry1: Float, _ s0: Int) -> (Float, Float, Float, Float) {
            let (v1x, v1y) = velocity(rx1, ry1, s0)
            let k1x = 2*h*v1x
            let k1y = 2*h*v1y
            
            let (v2x, v2y) = velocity(rx1+0.5*k1x, ry1+0.5*k1y, s2)
            let k2x = 2*h*v2x
            let k2y = 2*h*v2y
            
            let (v3x, v3y) = velocity(rx1+0.5*k2x, ry1+0.5*k2y, s2)
            let k3x = 2*h*v3x
            let k3y = 2*h*v3y
            
            let (v4x, v4y) = velocity(rx1+k3x, ry1+k3y, s4)
            let k4x = 2*h*v4x
            let k4y = 2*h*v4y
            
            return (v1x, v1y, (k1x+k4x+2.0*(k2x+k3x))/6.0, (k1y+k4y+2.0*(k2y+k3y))/6.0)
        }

        func smallRK1(_ rx1: Float, _ ry1: Float, _ v1x: Float, _ v1y: Float) -> (Float, Float) {
            /* Do the first small RK step.  No initial call to cart_velocity() is done
             * because it would be the same as the one above, so there's no need
             * to do it again */
            
            let k1x = h*v1x
            let k1y = h*v1y
            let (v2x, v2y) = velocity(rx1+0.5*k1x,ry1+0.5*k1y,s1)
            let k2x = h*v2x
            let k2y = h*v2y
            let (v3x, v3y) = velocity(rx1+0.5*k2x, ry1+0.5*k2y, s1)
            let k3x = h*v3x
            let k3y = h*v3y
            let (v4x, v4y) = velocity(rx1+k3x, ry1+k3y, s2)
            let k4x = h*v4x
            let k4y = h*v4y
            
            return ((k1x+k4x+2.0*(k2x+k3x))/6.0, (k1y+k4y+2.0*(k2y+k3y))/6.0)
        }
        
        func smallRK2(_ rx1: Float, _ ry1: Float, _ dx1: Float, _ dy1: Float) -> (Float, Float) {
        
            /* Do the second small RK step */
            
            let rx2 = rx1 + dx1
            let ry2 = ry1 + dy1
            
            let (v1x, v1y) = velocity(rx2, ry2, s2)
            let k1x = h*v1x
            let k1y = h*v1y
            let (v2x, v2y) = velocity(rx2+0.5*k1x, ry2+0.5*k1y, s3)
            let k2x = h*v2x
            let k2y = h*v2y
            let (v3x, v3y) = velocity(rx2+0.5*k2x, ry2+0.5*k2y, s3)
            let k3x = h*v3x
            let k3y = h*v3y
            let (v4x, v4y) = velocity(rx2+k3x, ry2+k3y, s4)
            let k4x = h*v4x
            let k4y = h*v4y
            
            return ((k1x+k4x+2.0*(k2x+k3x))/6.0, (k1y+k4y+2.0*(k2y+k3y))/6.0)
        }

        let s0 = s
        let s1 = (s+1) % SNAPSHOTCOUNT
        let s2 = (s+2) % SNAPSHOTCOUNT
        let s3 = (s+3) % SNAPSHOTCOUNT
        let s4 = (s+4) % SNAPSHOTCOUNT
        
        /* Calculate the density field for the four new time slices */
        
        density(t: t+0.5*h, s: s1)
        density(t: t+1.0*h, s: s2)
        density(t: t+1.5*h, s: s3)
        density(t: t+2.0*h, s: s4)
        
        /* Calculate the resulting velocity grids */
        
        vgrid(s1)
        vgrid(s2)
        vgrid(s3)
        vgrid(s4)
        
        /* Do all three RK steps for each point in turn */
        
        var esqmax: Float = 0.0
        var drsqmax: Float = 0.0
        
        let points2: [Point<Float>] = points.map { r in
                        
            let (v1x, v1y, dx12, dy12) = bigRKstep(r.x, r.y, s0)
            
            let (dx1, dy1) = smallRK1(r.x, r.y, v1x, v1y)
            
            let (dx2, dy2) = smallRK2(r.x, r.y, dx1, dy1)
            
            /* Calculate the (squared) error */
            
            let ex = (dx1+dx2-dx12)/15;
            let ey = (dy1+dy2-dy12)/15;
            esqmax = max(esqmax, ex*ex + ey*ey)
            
            /* Update the position of the vertex using the more accurate (two small
             * steps) result, and deal with the boundary conditions.  This code
             * does 5th-order "local extrapolation" (which just means taking
             * the estimate of the 5th-order term above and adding it to our
             * 4th-order result get a result accurate to the next highest order) */
            
            let dxtotal = dx1 + dx2 + ex   // Last term is local extrapolation
            let dytotal = dy1 + dy2 + ey   // Last term is local extrapolation
            drsqmax = max(drsqmax, dxtotal * dxtotal + dytotal * dytotal)
            
            return Point(x:clipf(f: r.x + dxtotal, min: 0, max: WIDTH),
                         y: clipf(f: r.y + dytotal, min: 0, max: HEIGHT))
        }
        
        return (sqrt(esqmax), sqrt(drsqmax), s4, points2)
    }
    
    
    /* Function to estimate the percentage completion */
    
    func complete(_ t: Float) -> Int {
        
        let res = 100.0 * logf(t/Float(INITH)) / logf(Float(EXPECTEDTIME/INITH))
        
        return res > 100.0 ? 100 : Int(res)
    }
    
    
    /* Function to do the transformation of the given set of points
     * to the cartogram */
    
    mutating func makeCartogram(blur: Float) -> [Point<Float>] {
        
        /* Create the grid of points */
        func create_grid() -> [Point<Float>] {
            var points = [Point<Float>]()
                            
            for row in 0...HEIGHT {
                for column in 0...WIDTH {
                    points.append(Point(x: Float(column), y: Float(row)))
                }
            }

            return points
        }

        /* Calculate the initial density and velocity for snapshot zero */
        
        density(t: 0.0, s: 0)
        vgrid(0)
        var s = 0
        
        /* Now integrate the points in the polygons */
        
        var step = 0
        var t = 0.5 * blur * blur
        var h = INITH
        
        var points = create_grid()
        
        var done: Int = 0
        repeat {
            
            /* Do a combined (triple) integration step */
            let error: Float
            let dr: Float
            let sp: Int
            (error, dr, sp, points) = twosteps(points: points, t: t, h: Float(h), s: s)
            
            /* Increase the time by 2h and rotate snapshots */
            
            t += Float(2.0 * h)
            step += 2
            s = sp
            
            /* Adjust the time-step.  Factor of 2 arises because the target for
             * the two-step process is twice the target for an individual step */
            
            let desiredratio = powf(TARGETERROR / error * 2.0, 0.2)
            if desiredratio > MAXRATIO {
                h *= MAXRATIO
            }
            else {
                h *= desiredratio
            }
            
            done = complete(t)
            print("Percent done: \(done)")
            
            /* If no point moved then we are finished */
            if dr <= 0.0 {
                break
            }
            
        } while true
        
        return points
    }
    
}
