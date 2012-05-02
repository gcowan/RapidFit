
#pragma once
#ifndef ObservableRef_H
#define ObservableRef_H

#include <vector>
#include <string>

using namespace::std;

class ObservableRef
{
	public:
		ObservableRef();
		ObservableRef( string );
		ObservableRef( vector<string> );
		~ObservableRef();

		operator string() const
		{
		        return Observable_Name;
		}
		operator vector<string>() const
		{
			return Observable_Names;
		}
		//operator const char *() const
		//{
		//	return Observable_Name.c_str();
		//}
		const ObservableRef operator []( unsigned int index )const
		{
			if( index < Observable_Refs.size() )
			{
				return Observable_Refs[index];
			}
			return ObservableRef("unknown");
		}

		void SetIndex( const int ) const;
		int GetIndex() const;
		string Name() const;
		string* NameRef() const;
		size_t size() const;
		void push_back( string );

	private:
		mutable string Observable_Name;
		mutable int Observable_Index;
		vector<string> Observable_Names;
		vector<ObservableRef> Observable_Refs;
};

#endif

