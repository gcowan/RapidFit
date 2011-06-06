#ifndef ObservableRef_H
#define ObservableRef_H

#include <vector>
#include <string>

using namespace::std;

class ObservableRef{

	public:
		ObservableRef();
		ObservableRef( string );
		ObservableRef( vector<string> );
		~ObservableRef();

		operator string() const
		{
		        string new_name = Observable_Name;
		        return new_name;
		}
		operator vector<string>() const
		{
			vector<string> new_names = Observable_Names;
			return new_names;
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

		void SetIndex( int );
		int GetIndex();
		string Name() const;
		string* NameRef();
		size_t size() const;
		void push_back( string );

	private:
		string Observable_Name;
		int Observable_Index;
		vector<string> Observable_Names;
		vector<ObservableRef> Observable_Refs;
};

#endif

