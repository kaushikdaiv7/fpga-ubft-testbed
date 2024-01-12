#include "host.h"
#include "constants.h"
#include <vector> 
#include <random>

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }
    
    clock_t htod, dtoh, comp; 

    /*====================================================CL===============================================================*/

    std::string binaryFile = argv[1];
    cl_int err;
    cl::Context context;
    cl::Kernel krnl1, krnl2;
    cl::CommandQueue q;
    
    auto devices = get_xil_devices();
    auto fileBuf = read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, 0, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            std::cout << "Setting CU(s) up..." << std::endl; 
            OCL_CHECK(err, krnl1 = cl::Kernel(program, "krnl", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    /*====================================================INIT INPUT/OUTPUT VECTORS===============================================================*/

    std::vector<uint32_t, aligned_allocator<uint32_t> > input(1);
    std::vector<uint8_t, aligned_allocator<uint8_t> > hash_hw(1);
    uint8_t *hash_sw = (uint8_t*) malloc(sizeof(uint8_t) * 1);

    

    /*====================================================SW VERIFICATION===============================================================*/

    /* Add software Blake*/

    /*====================================================Setting up kernel I/O===============================================================*/

    /* INPUT BUFFERS */
    OCL_CHECK(err, cl::Buffer buffer_input(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(data_t) * 1, input.data(), &err)); 

    /* OUTPUT BUFFERS */
    OCL_CHECK(err, cl::Buffer buffer_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(data_t) * 1, hash_hw.data(), &err)); 

    /* SETTING INPUT PARAMETERS */
    OCL_CHECK(err, err = krnl1.setArg(0, buffer_input));
    OCL_CHECK(err, err = krnl1.setArg(1, buffer_output));


    /*====================================================KERNEL===============================================================*/
    /* HOST -> DEVICE DATA TRANSFER*/
    std::cout << "HOST -> DEVICE" << std::endl; 
    htod = clock(); 
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_input}, 0 /* 0 means from host*/));
    q.finish();
    htod = clock() - htod; 
    
    /*STARTING KERNEL(S)*/
    std::cout << "STARTING KERNEL(S)" << std::endl; 
    comp = clock(); 
	OCL_CHECK(err, err = q.enqueueTask(krnl1));
    q.finish(); 
    comp = clock() - comp;
    std::cout << "KERNEL(S) FINISHED" << std::endl; 

    /*DEVICE -> HOST DATA TRANSFER*/
    std::cout << "HOST <- DEVICE" << std::endl; 
    dtoh = clock();
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();
    dtoh = clock() - dtoh;

    /*====================================================VERIFICATION & TIMING===============================================================*/

    printf("Host -> Device : %lf ms\n", 1000.0 * htod/CLOCKS_PER_SEC);
    printf("Device -> Host : %lf ms\n", 1000.0 * dtoh/CLOCKS_PER_SEC);
    printf("Computation : %lf ms\n",  1000.0 * comp/CLOCKS_PER_SEC);
    
    bool match = true;

    std::cout << "HW: " << hash_hw << " vs " << "SW: " << hash_sw << std::endl;
    if (hash_hw != hash_sw ){
        match = false; 
    }

    free(hash_sw);
    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << " count: " << count << std::endl;
    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}