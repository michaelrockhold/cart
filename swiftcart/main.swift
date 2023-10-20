import Foundation

let OFFSET: Float = 0.005

struct FilePopulation: Population {
    typealias ElementType = Float
    
    let height: Int
    let width: Int
    var data: [ElementType]
    
    init?() {
        guard let ln = readLine() else { return nil }
        var h = 1
        
        let tokens = ln.split(separator: " ")
        width = tokens.count
        data = tokens.map {
            return Float($0)!
        }
        while let ln = readLine() {
            h += 1
            let tokens = ln.split(separator: " ")
            data.append(contentsOf: tokens.map {
                return Float($0)!
            })
        }
        
        let sum = data.reduce(Double(0.0)){ (sum, elem) -> Double in
            return sum + Double(elem)
        }

        let mean = Float(sum / Double(data.count))
        
        data = data.map {
            $0 + OFFSET * mean
        }
        height = h
    }
}

func readpop() -> (Int, Int, [Float])? {
    guard let ln = readLine() else { return nil }
    var height = 1
    
    let tokens = ln.split(separator: " ")
    let width = tokens.count
    var data = tokens.map {
        return Float($0)!
    }
    while let ln = readLine() {
        height += 1
        let tokens = ln.split(separator: " ")
        data.append(contentsOf: tokens.map {
            return Float($0)!
        })
    }
    
    let sum = data.reduce(Double(0.0)){ (sum, elem) -> Double in
        return sum + Double(elem)
    }

    let mean = Float(sum / Double(data.count))
    
    data = data.map {
        $0 + OFFSET * mean
    }

    return (width, height, data)
}


func writepoints(gridx: [Float], gridy: [Float]) {
    for i in 0..<gridx.count {
        print("\(gridx[i]) \(gridy[i])")
    }
}


// Read in the population data
guard let rho = FilePopulation() else {
    print("Malformed input data")
    exit(-1)
}


// Allocate space for the cartogram code to use
// Transform the population data, then ignore it
// cart_makews(columnsCount, rowsCount);
//cart_transform(rho, columnsCount, rowsCount)
var cart = Cart.init(population: rho)

// Make the cartogram
//cart_makecart(gridx,gridy,(xsize+1)*(ysize+1),xsize,ysize,0.0);
let cgram = cart.makeCartogram(blur: 0.0)

// Write out the final positions of the grid points
for p in cgram {
    print("\(p.x) \(p.y)")
}
