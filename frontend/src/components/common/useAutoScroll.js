import { useEffect } from 'react';

const useAutoScroll = (dependency, ref) => {
    useEffect(() => {
        if (ref.current) {
            setTimeout(() => {
                ref.current.scrollIntoView({
                    behavior: 'smooth',
                    block: 'center'
                });
            }, 100);
        }
    }, [dependency, ref]);
};

export default useAutoScroll;
