/*
 * deformables.cpp
 *
 * Copyright (C) 2009 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "deformables.h"
#include "api/MegaMolCore.std.h"
#include "ModuleAutoDescription.h"
#include "vislib/vislibversion.h"
#include "vislib/Log.h"
#include "vislib/ThreadSafeStackTrace.h"

#include "SPHSimulation.h" 

/*
 * mmplgPluginAPIVersion
 */
DEFORMABLES_API int mmplgPluginAPIVersion(void) {
    return 100;
}


/*
 * mmplgPluginName
 */
DEFORMABLES_API const char * mmplgPluginName(void) {
    return "deformables";
}


/*
 * mmplgPluginDescription
 */
DEFORMABLES_API const char * mmplgPluginDescription(void) {
    return "This module helps simulate defomrabale objects using the SPH method";
}


/*
 * mmplgCoreCompatibilityValue
 */
DEFORMABLES_API const void * mmplgCoreCompatibilityValue(void) {
    static const mmplgCompatibilityValues compRev = {
        sizeof(mmplgCompatibilityValues),
        MEGAMOL_CORE_COMP_REV,
        VISLIB_VERSION_REVISION
    };
    return &compRev;
}


/*
 * mmplgModuleCount
 */
DEFORMABLES_API int mmplgModuleCount(void) {
    return 1; // TODO: Implement
}


/*
 * mmplgModuleDescription
 */
DEFORMABLES_API void* mmplgModuleDescription(int idx) {
	switch (idx) {
        case 0: return new megamol::core::ModuleAutoDescription<megamol::deformables::SPHSimulation>();
	}
    return NULL; // TODO: Implement
}


/*
 * mmplgCallCount
 */
DEFORMABLES_API int mmplgCallCount(void) {
    return 0; // TODO: Implement
}


/*
 * mmplgCallDescription
 */
DEFORMABLES_API void* mmplgCallDescription(int idx) {
    return NULL; // TODO: Implement
}


/*
 * mmplgConnectStatics
 */
DEFORMABLES_API bool mmplgConnectStatics(int which, void* value) {
    switch (which) {

        case 1: // vislib::log
            vislib::sys::Log::DefaultLog.SetLogFileName(static_cast<const char*>(NULL), false);
            vislib::sys::Log::DefaultLog.SetLevel(vislib::sys::Log::LEVEL_NONE);
            vislib::sys::Log::DefaultLog.SetEchoTarget(new vislib::sys::Log::RedirectTarget(static_cast<vislib::sys::Log*>(value)));
            vislib::sys::Log::DefaultLog.SetEchoLevel(vislib::sys::Log::LEVEL_ALL);
            vislib::sys::Log::DefaultLog.EchoOfflineMessages(true);
            return true;

        case 2: // vislib::stacktrace
            return vislib::sys::ThreadSafeStackTrace::Initialise(
                *static_cast<const vislib::SmartPtr<vislib::StackTrace>*>(value), true);

    }
    return false;
}
