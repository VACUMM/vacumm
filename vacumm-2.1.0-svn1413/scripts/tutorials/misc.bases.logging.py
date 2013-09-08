#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vacumm.misc.bases import Object

class MyObject(Object):
    pass

def main():
    
    MyObject.debug('Call logging classmethod debug')
    MyObject.info('Call logging classmethod info')
    MyObject.set_loglevel('debug')
    MyObject.get_logger().set_format('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    MyObject.debug('Call logging classmethod debug')
    
    # Show the default logging configuration
    obj = MyObject()
    obj.notice(
        'These are the default instance logging configuration:\n'
        '  level:       %r (%r)\n'
        '  format:      %r\n'
        '  date format: %r',
        obj.get_loglevel(), obj.get_logger().get_level(),
        obj.get_logger().get_format(),
        obj.get_logger().get_date_format()
    )
    
    # Configure logging at init
    obj = MyObject(
        logger_level='debug',
        logger_format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        logger_date_format='%H:%M:%S'
    )
    obj.notice(
        'Customizing instance logging configuration:\n'
        '  level:       %r (%r)\n'
        '  format:      %r\n'
        '  date format: %r',
        obj.get_loglevel(), obj.get_logger().get_level(),
        obj.get_logger().get_format(),
        obj.get_logger().get_date_format()
    )
    
    obj.debug('debug message')
    obj.verbose('verbose message')
    obj.info('info message')
    obj.notice('notice message')
    obj.warning('warning message')
    obj.error('error message')
    obj.critical('critical message')
    try: 0 / 0
    except Exception, e:
        obj.exception('Division by 0 failed !')
    
    MyObject.verbose('\n  MyObject.get_logger(): %r\n  MyObject.get_class_logger(): %r\n  MyObject.logger: %r', MyObject.get_logger(), MyObject.get_class_logger(), MyObject.logger)
    obj.verbose('\n  obj.get_logger(): %r\n  obj.get_class_logger(): %r\n  obj.logger: %r', obj.get_logger(), obj.get_class_logger(), obj.logger)
    
    obj.set_loglevel('notice')
    obj.notice('the loglevel is now %r', obj.get_loglevel())
    obj.verbose(
        'you will not see this message as it is emitted with a "verbose" level and '
        'the logger is now configured with a %r minimum level', obj.get_loglevel())
    
    obj.notice('Using the config method:\n  obj.logger.config(): %s\n  MyObject().logger.config(): %s\n  MyObject(logger_config=obj).logger.config(): %s',
        obj.logger.config(), MyObject().logger.config(), MyObject(logger_config=obj).logger.config())

if __name__ == '__main__':
    main()

