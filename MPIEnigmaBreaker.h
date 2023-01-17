
#include "EnigmaBreaker.h"

class MPIEnigmaBreaker : public EnigmaBreaker {
private:
    bool solutionFound( uint *rotorSettingsProposal );
    uint *expected;
    uint expectedLength;
public:
    MPIEnigmaBreaker( Enigma *enigma, MessageComparator *comparator );
    ~MPIEnigmaBreaker() override;
    void crackMessage() override;
    void getResult( uint *rotorPositions ) override;
    void setSampleToFind( uint *expected, uint expectedLength ) override;
    void numberToRotorPositions(uint64_t number, uint* rotorPositions);
};
