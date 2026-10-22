def exit_with_error(exit_code, log_message) {
    log.error log_message; exit exit_code
}


def checkForWhitespace(paramName, paramValue) {
    if (!paramValue) {
        return channel.empty()
    }

    return channel.fromPath(paramValue, glob: true)
        .map { path ->
            if (path.toString().contains(' ')) {
                // The 'error' command immediately terminates the pipeline with a message
                error "❌ The path from '${paramName}' contains whitespace: '${path}'. Please remove the spaces. ❌"
            }
            return path
        }
}