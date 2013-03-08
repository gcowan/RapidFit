/*!
 * @class ComponentRef
 * 
 * @brief This class has been written to store complex lookups performed when requesting a single component from a PDF
 * 
 * The important functions for a PDF developer are: get/setComponentNumber() and getComponentName()
 * 
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef _RAPIDFIT_COMPONENT_REF
#define _RAPIDFIT_COMPONENT_REF

///	System Headers
#include <string>

using namespace::std;

class ComponentRef
{
	public:
		/*!
		 * @brief Constructor with full name of component
		 *
		 * @param Name	This is the full name of the Component
		 *
		 * @param Obs	This is the name of the Observable being projected
		 */
		ComponentRef( const string Name, const string Obs );

		/*!
		 * @brief Destructor
		 */
		~ComponentRef();

		/*!
		 * @brief Correct Copy Constuctor for this class
		 */
		ComponentRef( const ComponentRef& );

		/*!
		 * @brief  Add a subComponetRef with the given name
		 *         i.e. add a new link in the chain
		 *
		 * @param Input   This is the Name of the new ComponentRef to be created and linked to this one
		 *
		 * @return Void
		 */
		void addSubComponent( const string );

		/*!
		 * @brief Get the subComponentRef for this instance
		 *        i.e. get the next link in the chain
		 *
		 * @warning This class will remove the object refrenced by this pointer on Destruction
		 *
		 * @return Returns the pointer to the CompoentRef internally stored by this object
		 */
		ComponentRef* getSubComponent() const;

		/*!
		 * @brief set the component number in this ComponentRef
		 *
		 * @oaram Input      This is a value you want to store after performing a lookup for this component in a list
		 *
		 * @return Void
		 */
		void setComponentNumber( const int Input );

		/*!
		 * @brief get the component number in this ComponentRef
		 *
		 * @return returns the stored result from a previous lookup
		 */
		int getComponentNumber() const;

		/*!
		 * @brief get the name of this ComponenRef
		 *
		 * @return returns the name of this component
		 */
		string getComponentName() const;

		string getObservableName() const;

		void setObservableName( const string );

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		ComponentRef& operator= ( const ComponentRef& );

		/*!
		 * This is the reference to the sub-ComponentRef stored in this instance
		 */
		ComponentRef* thisSubComponent;

		int thisIndex;		/*!	The index stored in this ComponentRef from an external lookup	*/
		string thisName;	/*!	The name of this particluar ComponentRef			*/
		string obsName;		/*!	The name of the Observable being Projected			*/
};

#endif

