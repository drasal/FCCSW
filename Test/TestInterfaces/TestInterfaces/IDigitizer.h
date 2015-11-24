//
//  IDigitizer.h
//  
//
//  Created by Julia Hrdinka on 14/04/15.
//
//

#ifndef IDIGITIZER_H
#define IDIGITIZER_H

#include "GaudiKernel/IAlgTool.h"

namespace fccedm {
  class ParticleCollection;
}

static const InterfaceID IID_IDigitizer ("IDigitizer", 1, 0);

class IDigitizer : virtual public IAlgTool {

public:
    static const InterfaceID& interfaceID() {
            return IID_IDigitizer;
    }
    
    virtual StatusCode smear(fccedm::ParticleCollection* particlecoll) = 0;
    
protected:
    ~IDigitizer() {}
};

#endif //IDIGITIZER_H
