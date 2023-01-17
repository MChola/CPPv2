
#include "math.h"
#include "MPIEnigmaBreaker.h"
#include "mpi/mpi.h"

MPIEnigmaBreaker::MPIEnigmaBreaker(Enigma *enigma, MessageComparator *comparator) : EnigmaBreaker(enigma, comparator) {
}

MPIEnigmaBreaker::~MPIEnigmaBreaker() {
    delete[] rotorPositions;
}

void MPIEnigmaBreaker::numberToRotorPositions(uint64_t number, uint *rotorPositions) {
    uint largestRotorSetting = enigma->getLargestRotorSetting();
    uint rotors = enigma->getNumberOfRotors();
    for (int i = 0; i < rotors; i++) {
        rotorPositions[i] = number % (largestRotorSetting + 1);
        number = number / (largestRotorSetting + 1);
    }
}

void informOtherProcesses (int size) {
    for (int j = 1; j < size; j++) {
        int doneFlag = 1;
        MPI_Request request;
        MPI_Isend(&doneFlag, 1, MPI_INTEGER, j, 400, MPI_COMM_WORLD, &request);
    }
}

void MPIEnigmaBreaker::setSampleToFind(uint *expected, uint expectedLength) {
    comparator->setExpectedFragment(expected, expectedLength);
    this->expectedLength = expectedLength;
    this->expected = expected;
}

bool MPIEnigmaBreaker::solutionFound(uint *rotorSettingsProposal) {
    for ( uint rotor = 0; rotor < rotors; rotor++ )
        rotorPositions[ rotor ] = rotorSettingsProposal[ rotor ];

    enigma->setRotorPositions(rotorPositions);
    uint *decodedMessage = new uint[ messageLength ];

    for (uint messagePosition = 0; messagePosition < messageLength; messagePosition++ )
        decodedMessage[ messagePosition ] = enigma->code(messageToDecode[ messagePosition ] );

    bool result = comparator->messageDecoded(decodedMessage);
    delete[] decodedMessage;
    return result;
}

void MPIEnigmaBreaker::getResult(uint *rotorPositions) {
    for ( uint rotor = 0; rotor < rotors; rotor++ ) {
        rotorPositions[ rotor ] = this->rotorPositions[ rotor ];
    }
}


void MPIEnigmaBreaker::crackMessage() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint rotorLargestSetting = enigma->getLargestRotorSetting();
    uint64_t numberOfCombinations = static_cast<uint64_t>(pow(rotorLargestSetting + 1, rotors));
    int numberOfTestProcessCombinations = (int) ((double) numberOfCombinations * 0.1) / size;

    MPI_Request mpi_r;
    MPI_Status mpi_s;
    int flag = 0;

    const uint BUFF_SIZE = sizeof(uint) * 999999;
    char buffer[BUFF_SIZE];
    int position = 0;

    if (!rank) {    ///////////////////////////////////////////////////////////////////////////////////////////////////0
        ///waiting for solution
        uint* goodPositions = new uint[rotors];
        MPI_Irecv(goodPositions, (int) rotors, MPI_UNSIGNED, MPI_ANY_SOURCE, 200, MPI_COMM_WORLD, &mpi_r);
        ///sending primary data
        MPI_Pack(&messageLength, 1, MPI_UNSIGNED, buffer, BUFF_SIZE, &position, MPI_COMM_WORLD);
        MPI_Pack(&expectedLength, 1, MPI_UNSIGNED, buffer, BUFF_SIZE, &position, MPI_COMM_WORLD);
        MPI_Pack(messageToDecode, (int) messageLength, MPI_UNSIGNED, buffer, BUFF_SIZE, &position, MPI_COMM_WORLD);
        MPI_Pack(expected, (int) expectedLength, MPI_UNSIGNED, buffer, BUFF_SIZE, &position, MPI_COMM_WORLD);
        MPI_Bcast(buffer, BUFF_SIZE, MPI_PACKED, MPI_ROOT_PROCESS_RANK, MPI_COMM_WORLD);

        ///time schedule
        MPI_Wtime();
        for (int i = 0; i < numberOfTestProcessCombinations; i++) {
            MPI_Test(&mpi_r, &flag, &mpi_s);
            if (flag) {
                rotorPositions = goodPositions;
                informOtherProcesses(size);
                return;
            }
            numberToRotorPositions(i, rotorPositions);
            if (solutionFound(rotorPositions)) {
                informOtherProcesses(size);
                return;
            }
        }
        double* times = new double[size];
        times[0] = MPI_Wtime();
        MPI_Status status;
        for (int i = 1; i < size; i++)
            MPI_Recv(&times[i], 1, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);

        uint combinationsCheckedDuringTest = numberOfTestProcessCombinations * size;
        double* combinationsPerSec = new double[size];
        double totalCombinationsPerSec = 0;
        for (int i = 0; i < size; i++) {
            combinationsPerSec[i] = numberOfTestProcessCombinations / times[i];
            totalCombinationsPerSec += combinationsPerSec[i];
        }
        double* speedPercentage = new double[size];
        for (int i = 0; i < size; i++)
            speedPercentage[i] = combinationsPerSec[i] / totalCombinationsPerSec;
        uint * combinationsScheduled = new uint [size];
        for (int i = 0; i < size; i++)
            combinationsScheduled[i] = (uint) (speedPercentage[i] * (double) (numberOfCombinations - combinationsCheckedDuringTest));

        uint* startNumber = new uint[size];
        uint* endNumber = new uint[size];   //end is start of next process
        for (int i = 0; i < size; i++) {
            if (i == 0)
                startNumber[i] = combinationsCheckedDuringTest;
            else
                startNumber[i] = endNumber[i-1];
            endNumber[i] = startNumber[i] + combinationsScheduled[i];
        }
        endNumber[size-1] = numberOfCombinations;    //last process gets all omitted numbers (due to rouding)
        for (int i = 1; i < size; i++) {
            MPI_Send(&startNumber[i], 1, MPI_UNSIGNED, i, 300, MPI_COMM_WORLD);
            MPI_Send(&endNumber[i], 1, MPI_UNSIGNED, i, 300, MPI_COMM_WORLD);
        }

        ///root process computing
        for (uint i = startNumber[0]; i < endNumber[0]; i++) {
            MPI_Test(&mpi_r, &flag, &mpi_s);
            if (flag) {
                rotorPositions = goodPositions;
                informOtherProcesses(size);
                return;
            }
            numberToRotorPositions(i, rotorPositions);
            if (solutionFound(rotorPositions)) {
                informOtherProcesses(size);
                return;
            }
        }
        MPI_Wait(&mpi_r, &mpi_s);
        rotorPositions = goodPositions;
    }

    else {  //////////////////////////////////////////////////////////////////////////////////////////////////////////!0
        ///waiting for solution
        MPI_Irecv(&flag, 1, MPI_INTEGER, MPI_ROOT_PROCESS_RANK, 400, MPI_COMM_WORLD, &mpi_r);
        ///receiving primary data
        MPI_Bcast(buffer, BUFF_SIZE, MPI_PACKED, MPI_ROOT_PROCESS_RANK, MPI_COMM_WORLD);
        MPI_Unpack(buffer, BUFF_SIZE, &position, &messageLength, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        MPI_Unpack(buffer, BUFF_SIZE, &position, &expectedLength, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        messageToDecode = new uint[messageLength];
        expected = new uint[expectedLength];
        MPI_Unpack(buffer, BUFF_SIZE, &position, messageToDecode, (int) messageLength, MPI_UNSIGNED, MPI_COMM_WORLD);
        MPI_Unpack(buffer, BUFF_SIZE, &position, expected, (int) expectedLength, MPI_UNSIGNED, MPI_COMM_WORLD);
        comparator->setMessageLength(messageLength);
        setSampleToFind(expected, expectedLength);

        ///time measure
        MPI_Wtime();
        for (int i = numberOfTestProcessCombinations * rank; i < numberOfTestProcessCombinations * (rank + 1); i++) {
            MPI_Test(&mpi_r, &flag, &mpi_s);
            if (flag)
                return;
            numberToRotorPositions(i, rotorPositions);
            if (solutionFound(rotorPositions)) {
                MPI_Isend(rotorPositions, (int) rotors, MPI_UNSIGNED, MPI_ROOT_PROCESS_RANK, 200, MPI_COMM_WORLD, &mpi_r);
                return;
            }
        }
        double time = MPI_Wtime();
        MPI_Send(&time, 1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);

        ///other processes computing
        uint startNumber, endNumber;
        MPI_Recv(&startNumber, 1, MPI_UNSIGNED, MPI_ROOT_PROCESS_RANK, 300, MPI_COMM_WORLD, &mpi_s);
        MPI_Recv(&endNumber, 1, MPI_UNSIGNED, MPI_ROOT_PROCESS_RANK, 300, MPI_COMM_WORLD, &mpi_s);
        for (uint i = startNumber; i < endNumber; i++) {
            MPI_Test(&mpi_r, &flag, &mpi_s);
            if (flag)
                return;
            numberToRotorPositions(i, rotorPositions);
            if (solutionFound(rotorPositions)) {
//                cout << rank << ": " << "FOUND!" << endl;
//                for (int j = 0; j < rotors; j++) {
//                    cout << rotorPositions[j] << endl;
//                }
                MPI_Isend(rotorPositions, (int) rotors, MPI_UNSIGNED, MPI_ROOT_PROCESS_RANK, 200, MPI_COMM_WORLD, &mpi_r);
                return;
            }
        }

    }

}
