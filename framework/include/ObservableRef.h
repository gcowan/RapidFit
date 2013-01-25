
#pragma once
#ifndef ObservableRef_H
#define ObservableRef_H

#include <vector>
#include <string>

using namespace::std;

class ObservableRef
{
	public:
		/*!
		 * @brief Constructor
		 */
		ObservableRef( string="UnNamed" );

		/*!
		 * @brief Copy Constructor
		 */
		ObservableRef( const ObservableRef& );

		/*!
		 * @brief Copy By Assignment
		 */
		ObservableRef& operator = ( const ObservableRef& );

		/*!
		 * @brief Destructor
		 */
		~ObservableRef();

		/*!
		 * @brief This Function casts the class into a string object
		 *
		 * @return Returns the name of the ObservableRef
		 */
		operator string() const
		{
		        return Observable_Name;
		}

		/*!
		 * @brief
		 *
		 * @param inputIndex
		 *
		 * @return void
		 */
		void SetIndex( const int ) const;

		/*!
		 * @brief
		 *
		 * @return 
		 */
		int GetIndex() const;

		/*!
		 * @brief
		 *
		 * @return
		 */
		string Name() const;

		/*!
		 * @brief
		 *
		 * @return
		 */
		const string* NameRef() const;

		/*!
		 * @brief Prints to screen the stored information internal to this class
		 *
		 * @return void
		 */
		void Print() const;

		/*!
		 * @brief This sets the unique ID based on the unique ID of some external object
		 *
		 * @return void
		 */
		void SetExternalID( size_t inputID ) const;

		/*!
		 * @brief This returns the unique ID of the class which is used to make the judgement as to whether the stored cache is applicable to the object calling it
		 *
		 * @return Returns a uniqueID which related to the system type as to how unique the object in memory can be
		 */
		size_t GetExternalID() const;

	private:
		string Observable_Name;
		mutable int Observable_Index;
		mutable size_t externID;
};

#endif

