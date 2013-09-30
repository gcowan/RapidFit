
#pragma once
#ifndef _THREADING_CONFIG_H_
#define _THREADING_CONFIG_H_

#include <string>

using namespace::std;

class ThreadingConfig
{
	public:
		string MultiThreadingInstance;
		unsigned int numThreads;
		ComponentRef* wantedComponent;
};

#endif

